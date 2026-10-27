# Mapping route: NHANES drugs → Decagon STITCH IDs

Research for issue #12 (map #1), point `4.2` of `PROJECT_CONTEXT.md`.
Date: 2026-10-27. Live checks were run against RxNav (`rxnav.nlm.nih.gov/REST`) on this date.

Status markers: **[verified]** = confirmed in a primary doc or a live API call; **[to verify]** = from memory or a secondary source, check before relying on it.

---

## Answer in one paragraph

Use the NHANES Multum code `RXDDRGID`, not the drug name. Prefix it with `GN` and look it up in RxNav (`findRxcuiById`, `idtype=MMSL_CODE`). This returns an RxNorm ingredient (`IN`) for single drugs and a multi-ingredient concept (`MIN`) for combination products. Split each `MIN` into its `IN`s with `related?tty=IN`. Then go from the `IN` to a PubChem CID (UNII from RxNorm → PubChem, with a name fallback), and normalise the CID to its PubChem parent (salt removed). Decagon's IDs look like `CID0xxxxxxxx`, which in STITCH's notation is the **stereo-specific / single-compound** form (`CIDs`), not the merged/flat form (`CID1…` = `CIDm`). So match on the exact PubChem CID first, and fall back to a **same-connectivity** match (stereo ignored) when the exact CID is missing. Log every fallback.

---

## 1. NHANES Multum code → RxNorm, without name matching

**NHANES side** [verified, CDC `RXQ_DRUG` documentation]
- `RXDDRGID` is the generic drug code. Each generic drug name has a unique code from the Cerner Multum Lexicon. Multum generic codes start with `d`.
- Codes starting with `a` are drugs NHANES participants reported that are **not in the Lexicon**. Codes starting with `c` are **unspecified products with a known therapeutic class**. Neither kind will map to RxNorm by code.
- The prescription file (`RXQ_RX_*`) links to the drug information file (`RXQ_DRUG`) on `RXDDRGID`.
- Source: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1988/DataFiles/RXQ_DRUG.htm

**RxNav side** [verified, RxNav `findRxcuiById` docs]
- `idtype=MMSL_CODE` is a *derived* identifier: the Multum MediSource Lexicon (SAB `MMSL`) term type plus the code. "The first two letters represent the term type: BD, BN, CD, GN or IN."
- Scope: `allsrc=0` (default) = Active concepts only (an RxNorm-sourced, unsuppressed atom). `allsrc=1` = Current concepts (any unsuppressed atom).
- Source: https://lhncbc.nlm.nih.gov/RxNav/APIs/api-RxNorm.findRxcuiById.html

**Live behaviour** [verified, 2026-10-27]

| Request | Result |
|---|---|
| `rxcui.json?idtype=MMSL_CODE&id=d00134` (bare NHANES code) | empty `idGroup`: the `GN` prefix is required |
| `rxcui.json?idtype=MMSL_CODE&id=GNd00134` | `6918` |
| `rxcui/6918/properties.json` | `metoprolol`, TTY `IN` |
| `rxcui/6918/allProperties.json?prop=all` | lists `MMSL_CODE = GNd00134` and `UNII_CODE = GEB06NHM23` (also DrugBank `DB00264`, ATC `C07AB02`) |
| `rxcui/214223/allProperties.json` (`amlodipine / benazepril`) | TTY `MIN`, `MMSL_CODE = GNd03829` |

**Consequences**
- Rule: `MMSL_CODE = "GN" + RXDDRGID`. One request per distinct code; NHANES has a few thousand distinct codes across all cycles, so the full lookup is small.
- The Multum generic code for metoprolol hangs on the salt-free `IN`, not on a salt `PIN`. The NHANES code is already at ingredient level, so salts mostly disappear at this step.
- Combination products have **their own** Multum `d` code (e.g. `d03829`), which maps to an RxNorm `MIN`. They are not stored as a list of single-ingredient codes.
- **Failure cases:** `a…` and `c…` codes; `d` codes that are old or retired in Multum and no longer carried in the current RxNorm MMSL source. Try `allsrc=1`, then fall back to name matching on `RXDDRUG` with `findRxcuiByString` / `approximateTerm`. Count each failure in a log. [The share of NHANES codes that fail is **to verify** on the real file.]
- For historical lookups, check `getRxcuiHistoryStatus` for concepts that were remapped. [to verify]

## 2. Splitting combinations; salts and esters

**RxNorm term types** [verified for the examples below; definitions per RxNorm TTY appendix, https://www.nlm.nih.gov/research/umls/rxnorm/docs/appendix5.html]
- `IN` = ingredient (base form, e.g. metoprolol). `PIN` = precise ingredient (salt / ester / specific form, e.g. metoprolol succinate). `MIN` = multiple ingredients (combination, e.g. amlodipine / benazepril).
- `MIN → IN`: `GET /REST/rxcui/{MIN}/related.json?tty=IN` returns the single ingredients (getRelatedByType, https://lhncbc.nlm.nih.gov/RxNav/APIs/api-RxNorm.getRelatedByType.html). [The live response for 214223 is **to verify**; a call failed for tooling reasons.]
- `PIN → IN`: `GET /REST/rxcui/{PIN}/related.json?tty=IN` returns the base ingredient (RxNorm relation `form_of` / `has_form`). Live: the name "metoprolol succinate" resolves to RxCUI `221124`. [That 221124 is a `PIN` with `form_of → 6918` is **to verify**.]
- On the MMSL-code route, NHANES codes land on `IN` or `MIN` directly (section 1), so the `PIN` step matters only for the name-matching fallback.
- **Edge cases:** some salts are separate `IN`s in RxNorm, not `PIN`s, when the salt is clinically distinct (e.g. some potassium / sodium salts, or esters used as distinct drugs). Also, if a Decagon drug is itself registered at the salt CID in STITCH, stripping the salt would lose the match. Section 4 normalises both sides to the same level.

## 3. RxNorm ingredient → PubChem CID

**Options**
1. **UNII (preferred).** RxNorm `IN` concepts carry a `UNII_CODE` property [verified: metoprolol → `GEB06NHM23`]. FDA's GSRS deposits substances into PubChem keyed by UNII, and PubChem compound synonyms include the UNII. Lookup: `GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{UNII}/cids/JSON`. [The PubChem response is **to verify** live.] A UNII is substance-level and stereo-specific, so the CID it returns is specific (e.g. a racemate vs. a single enantiomer).
2. **Name (fallback).** `compound/name/{IN name}/cids/JSON` uses PubChem's synonym index. Common names can hit several CIDs or the wrong form; take the first CID only with a log entry.
3. **DrugBank ID (cross-check).** RxNorm also carries DrugBank IDs (e.g. `DB00264`), useful for checking a mapping; DrugBank's own data needs an academic licence.

**Normalisation in PubChem PUG-REST** (https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) [endpoint names **to verify**]
- Parent compound (salt or solvate removed): `compound/cid/{cid}/cids/JSON?cids_type=parent`.
- Same skeleton with stereo ignored: `compound/fastidentity/cid/{cid}/cids/JSON?identity_type=same_connectivity`, or compare the first 14-character block of the InChIKey (`compound/cid/{cid}/property/InChIKey/JSON`).

**PubChem rate limits** (https://pubchem.ncbi.nlm.nih.gov/docs/programmatic-access) [**to verify** exact wording]: at most 5 requests per second, 400 per minute, and 300 s of running time per minute. Throttling is dynamic and shown in response headers. With a few hundred ingredients this is not a bottleneck; cache everything.

## 4. Are Decagon's IDs flat (`CIDm`) or stereo (`CIDs`)?

- STITCH has two kinds of chemical ID. Older files (STITCH ≤ 4 and the derived SIDER/TWOSIDES/Decagon files) use `CID0…` for **single / stereo-specific** compounds and `CID1…` for **merged (flat)** compounds. STITCH 5 renamed these `CIDs…` and `CIDm…`. In both, the digits after the flag are a PubChem CID. Source: a reply in Biostars thread https://biostars.org/p/155346 (about 2015), which says "CID0s are individual compounds, CID1s are the merged compounds" and announces the rename to `CIDs` / `CIDm`. [Stated as coming from the STITCH team; checking against the STITCH 5 paper, Szklarczyk et al. 2016 NAR, and the STITCH download README is **to verify**. stitch.embl.de could not be reached during this research because of a TLS certificate error.]
- Decagon IDs have the shape `CID` + 9 digits, e.g. `CID000002173` = `CID` + flag `0` + `00002173`. Under the convention above that is the **stereo / single-compound** form, with PubChem CID `2173`. [That *all* 645 Decagon IDs start with `CID0` should be checked on the downloaded file: `cut -d, -f1,2 bio-decagon-combo.csv | grep -c CID1`.]
- In SIDER the same drug appears as both, with **different numbers**: the flat ID points to a representative CID with stereo removed, the stereo ID to the specific compound. So the flag is not cosmetic; the CID numbers themselves differ. [to verify on SIDER files]
- Decagon strips nothing itself; its drug IDs come from TWOSIDES (Tatonetti et al. 2012). The paper does not say flat vs. stereo. [Check the Decagon paper's data section.]

**What it means for matching**
- Exact match: the PubChem CID from section 3 (after the parent step) equals the Decagon number. This works when STITCH's chosen compound and the UNII-linked compound are the same stereo form.
- Mismatch risk: drugs sold as racemates vs. a single enantiomer (e.g. esomeprazole / omeprazole, levocetirizine / cetirizine, escitalopram / citalopram), and drugs whose PubChem record exists with and without defined stereocentres. In these cases the exact CIDs differ.
- Recommended fallback: map *both* sides to the same-connectivity level (InChIKey first block, or `fastidentity … same_connectivity`) and match there. Precompute this for Decagon's 645 CIDs once. **Caution:** at the connectivity level, enantiomer pairs (omeprazole / esomeprazole) collapse to the same key. Report such many-to-one matches rather than merging them silently.

## Recommended route (summary)

1. Unique `RXDDRGID` across the pooled NHANES cycles → drop `a…` / `c…` codes into the "unmappable by code" log.
2. `GN + RXDDRGID` → RxNav `findRxcuiById` (`MMSL_CODE`) → RxCUI. On failure: `allsrc=1`, then name match on `RXDDRUG`.
3. If TTY = `MIN` → `related?tty=IN`. If `PIN` → `related?tty=IN`. Result: one or more `IN` RxCUIs.
4. `IN` → `UNII_CODE` (RxNav property) → PubChem CID. Fallback: `IN` name → PubChem CID.
5. CID → PubChem parent CID.
6. Match against Decagon's 645 numbers (strip `CID` and the leading flag digit / zeros): exact first, then same-connectivity. Record the match level per drug.
7. Cache every response (RxNav asks for 12–24 h caching; for a fixed analysis, freeze to disk and record the RxNorm version from `/REST/version.json`).

**Rate limits**
- RxNav: "no more than 20 requests per second per IP address", caching for 12–24 h recommended; 100 req/s with a UTS API key; RxNav-in-a-Box (Docker) for bulk use. [verified, https://lhncbc.nlm.nih.gov/RxNav/TermsofService.html]
- PubChem PUG-REST: 5 req/s, 400 req/min. [to verify]

**Known failure cases (to log and count)**
- NHANES `a…` (not in Multum) and `c…` (unspecified, class only) codes.
- Retired Multum codes not present in the current RxNorm MMSL source.
- Combination products where one component is not a chemical (vitamins, multi-electrolyte solutions, herbal / biologic ingredients): no meaningful PubChem CID.
- Biologics (insulins, monoclonal antibodies): no small-molecule CID; Decagon (STITCH) has none of these anyway.
- Racemate vs. enantiomer CID mismatches (section 4).
- Salts registered as separate `IN`s, and Decagon drugs registered at a salt CID.
- UNII missing in RxNorm, or UNII not yet linked to a PubChem compound (a substance without a standardised structure).
