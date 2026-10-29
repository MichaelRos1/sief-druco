# Decagon files: columns, counts and coverage

Facts for issue #4 (map #1). Checked 2026-10-29 on the downloaded files. Point IDs (`3.1`, `3.3.1`, ...) refer to `PROJECT_CONTEXT.md`. Nothing is decided here; the facts feed #5 (graph inputs) and #6 (splits).

Every number below was computed from the files with pandas (all columns read as strings). Nothing here comes from memory.

## 1. Short answer

| Question from #4 | Answer | `PROJECT_CONTEXT.md` said |
|---|---|---|
| Column names | As listed in `3.1`, exactly (section 3) | same, marked "verify" |
| Side effects with ≥ 500 pairs | **963**, not 964. No counting rule gives 964 on this file (section 4) | 964 |
| Combo drugs with target edges, `targets.csv` | **284** of 645 (361 without) | 284 / 361 |
| Same, `targets-all.csv` | **284** of 645. It adds **no** drug and **no** edge for the 645 (section 6) | open question (`3.3.4`) |
| Combo drugs with a mono profile | **639** of 645 (6 without) | 640 |
| Target proteins missing from the PPI graph | **8** proteins (94 edges) in `targets.csv`. No drug loses all its targets | open |

## 2. Download

- **Source:** `https://snap.stanford.edu/decagon/bio-decagon-<name>.tar.gz` for `combo`, `mono`, `targets`, `targets-all`, `ppi`, `effectcategories`.
- **Date:** `combo` downloaded 2026-10-28, the other five on 2026-10-29. Each archive holds one CSV. The archive dates are 2018-02-04, except `targets-all`, which is dated 2018-02-16.
- **Location:** `local/data/raw/decagon/` (git-ignored), with `SHA256SUMS` next to the files.

| File | SHA-256 of the `.tar.gz` | SHA-256 of the `.csv` |
|---|---|---|
| combo | `ab5a774019fb39da6434fdde0a3827bc02f7b54fbcedfabd945450ebd4030f59` | `60bd3d071345473236f56638ebb732c3c964f0fd3468fdd29a444bcc1bf80d35` |
| mono | `26e2200092d4e974f1cfae1ce23446ce42ff31a13630b6e0a170153d57451be7` | `c82583ab8159ba42f8c771adfbc403f4243ec502dc85f732b0f62023b20aeb74` |
| targets | `d3c561a00edf1a174ebbf2d7b5c8fe4e2fcfb249f31d3b747de46da645613d68` | `dd7ca7ecbe4026c6be15c043d9214ec30ad36db32b2ccceeae4091b3b9f1e564` |
| targets-all | `a9cee9710722802287e76a63cbd5647686d4f93f3f43ab534ce2bf264f9086b8` | `14f26ade2e1376acdf38960ed6568da2b41bf33dde114dbcde69218738a09110` |
| ppi | `d9b5e5e42c38aaacf5ba2c55ee8db27551626b9cd1296757b29c932d916519fc` | `3665bf9b78b1d0a99c2504879eafad06cb4cf708aecbcac0c480e90e5e77027c` |
| effectcategories | `1fc7d66aa0fcac404b707a661c20b790200a2c7ff6b404f40b11b8e0e6840e38` | `7841f4df9441e903618d42f520cb8196adf6aad6b9cf4d5b016c684500f69437` |

## 3. Columns and row counts

No file has duplicate rows or empty cells.

| File | Columns | Rows | Unique IDs |
|---|---|---|---|
| `bio-decagon-combo.csv` | `STITCH 1, STITCH 2, Polypharmacy Side Effect, Side Effect Name` | 4,649,441 | 645 drugs; 63,473 pairs; 1,317 side effects |
| `bio-decagon-mono.csv` | `STITCH, Individual Side Effect, Side Effect Name` | 174,977 | 639 drugs; 10,184 side effects |
| `bio-decagon-targets.csv` | `STITCH, Gene` | 18,690 | 284 drugs; 3,648 proteins |
| `bio-decagon-targets-all.csv` | `STITCH, Gene` | 131,034 | 1,774 drugs; 7,795 proteins |
| `bio-decagon-ppi.csv` | `Gene 1, Gene 2` | 715,612 | 19,081 proteins |
| `bio-decagon-effectcategories.csv` | `Side Effect, Side Effect Name, Disease Class` | 561 | 561 side effects; 37 classes |

Notes on structure that the data layer (`5.1`) relies on:

- **Combo pairs are listed once, in one order.** No pair appears as both (a, b) and (b, a), no row has a = b, and no (pair, side effect) repeats. Ordered and unordered pair counts are both 63,473. Reverse edges must be added in code.
- **PPI edges are listed once, in one order.** 715,612 undirected edges, no self loops, no row in both directions. Reverse edges must be added in code.
- **Side-effect names are 1 : 1 with IDs** in `combo` (1,317 IDs, 1,317 names).

## 4. The ≥ 500-pair filter: 963, not 964

Rule: per `Polypharmacy Side Effect`, count distinct drug pairs; keep side effects with ≥ 500.

- **Result: 963 side effects.** Thresholds 499, 500, 501 and 502 all give 963. The sorted counts near the cut are …, 492, 492, **498**, **502**, 504, …. Only a threshold of 493–498 gives 964, which would add avascular necrosis (`C0221247`, 498 pairs).
- **No other counting rule changes this.** The authors' loader counts ordered (drug 1, drug 2) keys per side effect (`research/severity-weights.md` §6). On this file ordered keys, unordered pairs and raw rows are all the same number, because no pair is listed twice (section 3). So the earlier open point in `severity-weights.md` §6 ("whether counting raw rows instead … gives 964") is now settled: it does not.
- **Likely cause (not verified):** the paper describes 1,318 side-effect types and 4,651,131 associations. The released file has 1,317 and 4,649,441, that is 1 type and 1,690 rows fewer. The paper's 964 was probably computed on that slightly larger version. That version is not public, so the gap can't be closed.
- **After the filter:** 4,576,287 triples (98.4%), still 645 drugs, **63,472 pairs**. One pair, `CID000003914 | CID000005379`, loses its only label (iridocyclitis) and drops out of the labelled set.
- **Labels per pair after the filter:** mean 72.10, median 52, SD 64.21, min 1, max 505. The `3.3.3` figures (73.25 / 53 / 66.05) are for all 1,317 types; they match the BioSNAP page, which is computed on the unfiltered file.
- **Density after the filter:** 63,472 of 207,690 possible pairs labelled (30.56%). Ones within labelled pairs: 7.49%.
- **963 has precedent:** Malone et al. 2018 report 963 side effects (and 4,649,441 interactions, matching this file), and Lloyd et al. 2024 report a median over 963 (`research/decagon-stricter-splits.md`).
- Largest side effects (names as in the file): "arterial pressure NOS decreased" (`C0020649`, 28,568 pairs), "anaemia" (`C0002871`, 27,006), "Difficulty breathing" (`C0013404`, 26,037). Smallest kept: 502 pairs.

## 5. Mono profiles

- **639 drugs, 10,184 side-effect IDs** (`3.1` says 640 / 10,185). Every mono drug is a combo drug.
- **6 combo drugs have no mono profile:** `CID000002673`, `CID000003562`, `CID000003928`, `CID000004723`, `CID000005052`, `CID000123620`. With the planned input layer they get an all-zero feature vector.
- Side effects per drug: mean 273.8, median 230, min 5, max 1,550.
- **The mono vector is very sparse across drugs.** Of the 10,184 mono side effects, 2,032 occur in only one drug. Counts present in at least k drugs:

  | k | 1 | 2 | 5 | 10 | 20 | 50 |
  |---|---|---|---|---|---|---|
  | mono side effects | 10,184 | 8,152 | 5,640 | 3,858 | 2,373 | 934 |

- **Mono and combo side-effect IDs don't overlap at all.** None of the 1,317 combo CUIs appears as a mono CUI. Even where names match (184 of 1,317 combo names, case-insensitive), the CUIs differ. For example, combo "Blood calcium decreased" is `C0020598`, while mono has the same name under four other CUIs. So mono features can't be matched to pair labels by ID. This matters for the per-drug additive baseline (`6.5.1`), which therefore has to learn its own per-drug, per-side-effect term rather than copy a mono flag.

## 6. Targets

- **`targets.csv`:** 18,690 edges, all on combo drugs. **284** of the 645 have at least one target; **361** have none. Targets per covered drug: mean 65.8, median 11, max 2,840.
- **`targets-all.csv`:** 131,034 edges on 1,774 drugs, but 1,490 of those drugs are outside the 645. Restricted to the 645, it holds **exactly the same 18,690 edges** on the same 284 drugs. Every `targets.csv` edge is in `targets-all.csv`. **For this project the two files are the same**; `targets-all` gives no extra coverage of the 361 (`3.3.4` answered).
- **Proteins missing from the PPI graph:** 8 of `targets.csv`'s 3,648 proteins (94 edges) are not among the 19,081 PPI proteins. Every covered drug still has at least one target in the PPI graph. In `targets-all.csv`, 41 proteins are missing, but none of the extra ones touch the 645.
- **Drugs with no signal from either source:** 5 combo drugs have neither a target nor a mono profile (5 of the 6 drugs in section 5). They enter the model only through their drug–drug edges.

## 7. Effect categories

`effectcategories.csv` covers **399 of the 963** filtered side effects (561 rows, 37 disease classes). This agrees with `research/severity-weights.md`.

## 8. What this means for later tickets (facts only)

- **#5 Graph inputs:** the choice of target file is moot: both files give the same 18,690 edges on 284 drugs. Still open there: how to handle the 361 drugs with no targets, the 6 with no mono profile and the 5 with neither, and whether to filter the 10,184-dim mono vector (2,032 columns are single-drug).
- **#5 / `5.1`:** the side-effect set is 963 under the stated rule. The data layer either uses 963 (and RGCN uses `num_relations = 966`) or explicitly reproduces 964 by another rule. Pick one and say so when comparing with published numbers.
- **#6 Splits:** 63,472 labelled pairs after the filter, not 63,473. Pairs and PPI edges are stored one way only, so the split must be done on unordered pairs and reverse edges added afterwards.
- **`3.2` node count:** 19,081 PPI proteins, plus 8 target-only proteins if they are kept as nodes (19,089). Malone et al. 2018 report 19,089 proteins, consistent with keeping them (`research/decagon-stricter-splits.md`).
