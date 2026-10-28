# Mapping route: NHANES drugs → Decagon STITCH IDs

Research for issue #12 (map #1), point `4.2` of `PROJECT_CONTEXT.md`.
Date: 2026-10-27, updated 2026-10-28 after the route decision on the issue ("Route decided: reverse direction"). Live checks were run against RxNav (`rxnav.nlm.nih.gov/REST`, data version `05-Oct-2026`, API `3.1.356` per `/REST/version.json`) and PubChem PUG-REST / PUG-View on these dates.

Status markers: **[verified]** = confirmed in a primary doc or a live API call; **[to verify]** = from memory or a secondary source, check before relying on it; **[corrected]** = a preliminary claim from 2026-10-27 that turned out wrong; **[blocked]** = could not be checked, see *Blockers* at the end.

---

## Answer in one paragraph

Meet the two sides at the RxNorm ingredient (`IN`) level. **NHANES side:** prefix the Multum code `RXDDRGID` with `GN` and look it up in RxNav (`findRxcuiById`, `idtype=MMSL_CODE`); split `MIN`s and lift `PIN`s with `related.json?tty=IN` (both verified live). **Decagon side (645 drugs, once):** strip `CID` and the flag digit to get a PubChem CID, find a UNII for it in PubChem, and look the UNII up in RxNav with `findRxcuiById?idtype=UNII_CODE` (verified live). Decagon-style CIDs are often salt or stereo-undefined records with no UNII of their own (e.g. CID 2173, a stereo-undefined ampicillin anion), so the hop needs a fallback chain: exact CID → PubChem parent CID → same-connectivity CIDs (batched synonyms, `UNII-…`) → name match in RxNav. The STITCH README confirms `CID0…` / `CIDs` = "stereo-specific" and `CID1…` / `CIDm` = "flat", but the one Decagon ID checked (`CID000002173`) carries the number that SIDER uses as the **flat** ID of ampicillin, so the earlier conclusion "Decagon IDs are stereo-specific" does not hold up [corrected]. Whether all 645 Decagon IDs start with `CID0` is still open: the combo file is 35.7 MB [blocked].

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
- **Failure cases:** `a…` and `c…` codes; `d` codes that are old or retired in Multum and no longer carried in the current RxNorm MMSL source. Try `allsrc=1`, then fall back to name matching on `RXDDRUG` with `findRxcuiByString` / `approximateTerm`. Count each failure in a log. [The share of NHANES codes that fail is **blocked**: the prescription file `RXQ_RX_J.xpt` is 9.6 MB, see *Blockers*.]

**Retired concepts: `getRxcuiHistoryStatus`** [verified, docs https://lhncbc.nlm.nih.gov/RxNav/APIs/api-RxNorm.getRxcuiHistoryStatus.html and live calls 2026-10-28]
- Request: `GET /REST/rxcui/{rxcui}/historystatus.json`. Scope "Current and Historical". Status values per the docs: `Active`, `Obsolete`, `Remapped`, `Quantified`, `NotCurrent`, `Unknown` (never in any monthly release).
- For `Remapped` (and `NotCurrent` with a designated successor) the replacement is in `derivedConcepts.remappedConcept[]` with `remappedRxCui`, `remappedName`, `remappedTTY`, `remappedActive`. `derivedConcepts.ingredientConcept[]` gives the `IN`s of a concept.
- The docs point to `findActiveProducts` as the more convenient call for translating an inactive RxCUI to active concepts.
- Live:

| Request | Result |
|---|---|
| `rxcui/6918/historystatus.json` | `status: Active`, `isCurrent: YES`, `activeStartDate: 042005`, TTY `IN` |
| `rxcui/1306116/historystatus.json` (`Oxiglutatione`, an `IN` taken from `/REST/allstatus.json?status=Remapped`) | `status: Remapped`, `remappedDate: 122013`, `remappedConcept` → `25953 glutathione disulfide (IN, active YES)` |
| `rxcui/1309219/historystatus.json` (`poppyseed oil`, `IN`) | `Remapped` → `1432299 poppy seed oil (IN)` |

- Use on the NHANES side: if an MMSL lookup returns an RxCUI that is not `Active`, follow `remappedConcept` to the current `IN`. Note the input here is an RxCUI; a Multum code that is absent from current RxNorm gives no RxCUI at all, so history status does not help there (name fallback does).

## 2. Splitting combinations; salts and esters

**RxNorm term types** [verified for the examples below; definitions per RxNorm TTY appendix, https://www.nlm.nih.gov/research/umls/rxnorm/docs/appendix5.html]
- `IN` = ingredient (base form, e.g. metoprolol). `PIN` = precise ingredient (salt / ester / specific form, e.g. metoprolol succinate). `MIN` = multiple ingredients (combination, e.g. amlodipine / benazepril).
- `MIN → IN` [verified live 2026-10-28]: `GET /REST/rxcui/214223/related.json?tty=IN` → `17767 amlodipine (IN)` and `18867 benazepril (IN)`. `rxcui/214223/properties.json` confirms TTY `MIN`. (getRelatedByType, https://lhncbc.nlm.nih.gov/RxNav/APIs/api-RxNorm.getRelatedByType.html)
- `PIN → IN` [verified live 2026-10-28]: `rxcui/221124/properties.json` → `metoprolol succinate`, TTY `PIN`. `rxcui/221124/related.json?tty=IN` → `6918 metoprolol (IN)`. `rxcui/221124/related.json?rela=form_of` → also `6918 metoprolol (IN)`, so the link is the `form_of` relation.
- On the MMSL-code route, NHANES codes land on `IN` or `MIN` directly (section 1), so the `PIN` step matters only for the name-matching fallback. **On the Decagon side the `PIN` step is needed too:** a salt CID can resolve to a `PIN` (live: CID 60496 amlodipine besylate → UNII `864V2Q084H` → `104416 amlodipine besylate (PIN)`), so every RxCUI reached from PubChem goes through `related.json?tty=IN` as well.
- **Edge cases:** some salts are separate `IN`s in RxNorm, not `PIN`s, when the salt is clinically distinct (e.g. some potassium / sodium salts, or esters used as distinct drugs). [to verify: no example checked live.]

## 3. Reverse hop: PubChem CID → RxNorm ingredient (Decagon side)

This replaces the forward direction (RxNorm `IN` → UNII → PubChem CID) from the first version of these notes. The forward lookup `compound/name/{UNII}/cids` was never tested live and is no longer needed.

**Which RxNav endpoint accepts UNII** [verified]
- `GET /REST/idtypes.json` lists `UNII_CODE` among the id types (with `AMPID`, `ATC`, `DRUGBANK`, `MMSL_CODE`, `NDC`, `SNOMEDCT`, …).
- The `findRxcuiById` docs list `UNII_CODE` as "FDA Unique Ingredient Identifier Code" (example `LVX8N1UT73`). https://lhncbc.nlm.nih.gov/RxNav/APIs/api-RxNorm.findRxcuiById.html
- Live: `GET /REST/rxcui.json?idtype=UNII_CODE&id=GEB06NHM23` → `6918` (metoprolol), same with `allsrc=1`.

**CID → UNII in PubChem** [verified live, two ways]
- PUG-View, one CID per request: `GET https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON?heading=UNII`. The value sits under *Names and Identifiers → Other Identifiers → UNII*, sourced from "FDA Global Substance Registration System (GSRS)". Returns HTTP 404 when the compound has no UNII.
- PUG-REST synonyms, batched: `POST https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/synonyms/JSON` with `cid=<comma list>`. UNIIs appear as synonyms of the form `UNII-<code>` (live: CID 6249 → `UNII-7C782967RD`, CID 4171 → `UNII-GEB06NHM23`). This is the cheaper way for many CIDs.

**Live test of the full hop** (2026-10-28; CIDs looked up by `compound/name/{name}/cids/JSON`: metoprolol 4171, amlodipine 2162, aspirin 2244, warfarin 54678486)

| PubChem CID | PubChem title (InChIKey) | UNII | RxNav `UNII_CODE` → RxCUI |
|---|---|---|---|
| 4171 | Metoprolol | `GEB06NHM23` | `6918 metoprolol (IN)` |
| 2162 | Amlodipine | `1J444QC288` | `17767 amlodipine (IN)` |
| 2244 | Aspirin | `R16CO5Y76E` | `1191 aspirin (IN)` |
| 54678486 | Warfarin | `5Q7ZVV76EI` | `11289 warfarin (IN)` |
| 2771 | Citalopram | `0DHU5B8D6V` | `2556 citalopram (IN)` |
| 4594 | Omeprazole | `KG60484QX9` | `7646 omeprazole (IN)` |
| 5291 | Imatinib | `BKJ8M8G5HI` | `282388 imatinib (IN)` |
| 60496 | Amlodipine Besylate | `864V2Q084H` | `104416 amlodipine besylate (PIN)` → needs `PIN → IN` |
| **2173** (from `CID000002173`) | IUPAC name of the ampicillin anion, `AVKUERGKIZMTKX-UHFFFAOYSA-M` (no stereo, charged) | **none** (PUG-View 404, no `UNII-` synonym) | direct hop fails |
| 5001396 (SIDER flat ID of warfarin) | warfarin sodium, `PJVWKTKQMONHTI-UHFFFAOYSA-M` | none | direct hop fails |

**Fallback chain** [verified live on CID 2173 and 5001396]
1. Exact CID → UNII → RxNav.
2. Parent CID: `GET /rest/pug/compound/cid/{cid}/cids/TXT?cids_type=parent` → UNII → RxNav. Live: 5001396 → parent `54678486` → `5Q7ZVV76EI` → `11289 warfarin (IN)`. For 2173 the parent is `2174` (neutral ampicillin, still stereo-undefined), which has no UNII either.
3. Same connectivity (stereo ignored): `GET /rest/pug/compound/cid/{parent}/cids/TXT?cids_type=same_connectivity` (or `compound/fastidentity/cid/{cid}/cids/TXT?identity_type=same_connectivity`, same 47 CIDs for 2174), then batched synonyms for all of them. Live for 2174: 2 of 47 CIDs carry a UNII: `6249 ampicillin` → `7C782967RD` → `733 ampicillin (IN)`; `448911 L-(+)-Ampicillin` → `SFV84J2559` → no RxNorm concept. Keep the hits that reach RxNorm; if more than one distinct `IN` is reached, log it for a hand check (enantiomer pairs land here).
4. Name match: PubChem title → RxNav `GET /REST/rxcui.json?name={name}&search=2` (findRxcuiByString, normalized search) or `GET /REST/approximateTerm.json?term={name}&maxEntries=3&option=1`. Live: the title is useless when PubChem's title is an IUPAC name (2173, 2174 gave no hit); for normal titles it agrees with the UNII route on all 7 drugs above. `approximateTerm` on `Dactinomycin (USP)` → rank-1 `3100`, which `findRxcuiByString name=dactinomycin` also returns. Docs: https://lhncbc.nlm.nih.gov/RxNav/APIs/api-RxNorm.getApproximateMatch.html
5. Every RxCUI reached → `related.json?tty=IN` and `historystatus.json` if not active. [verified live: on an `IN` it returns the `IN` itself (`6918/related.json?tty=IN` → `6918 metoprolol`); on the `PIN` from CID 60496 it lifts to the base (`104416/related.json?tty=IN` → `17767 amlodipine`)]. So one call works for `IN`, `PIN` and `MIN` alike.

**PubChem endpoint names** [verified in the PUG-REST docs, Markdown version https://pubchem.ncbi.nlm.nih.gov/pcfe/docs/markdown/pug-rest.md, and live]
- `cids_type` accepts `original, parent, component, preferred, same_stereo, same_isotopes, same_connectivity, same_tautomer, same_parent, same_parent_stereo, same_parent_isotopes, same_parent_connectivity, same_parent_tautomer`.
- `fastidentity` takes `identity_type` of `same_connectivity, same_tautomer, same_stereo, same_isotope, same_stereo_isotope, nonconflicting_stereo, same_isotope_nonconflicting_stereo`.
- The preliminary form `compound/cid/{cid}/cids/JSON?cids_type=parent` was right; `fastidentity/cid/{cid}/cids/JSON?identity_type=same_connectivity` also works.

**PubChem rate limits** [verified, https://pubchem.ncbi.nlm.nih.gov/pcfe/docs/markdown/programmatic-access.md ("Last updated: May 19, 2025") and https://pubchem.ncbi.nlm.nih.gov/pcfe/docs/markdown/dynamic-request-throttling.md ("Last updated: August 19, 2026")]
- Usage policy: no more than **5 requests per second**; no API keys or whitelists to exceed this.
- Throttling is dynamic; requests over the limit get **HTTP 503**, and repeat offenders are blocked for a period of time. Every PUG-REST response carries `X-Throttling-Control` (Request Count / Request Time / Service status, Green < 50 %, Yellow 50–75 %, Red > 75 %, Black = blocked). Live example: `Request Count status: Green (1%), Request Time status: Green (0%), Service status: Green (6%)`.
- **[corrected]** The preliminary "400 per minute, 300 s of running time per minute" is not in the current docs (neither page mentions a per-minute figure). Only 5 req/s is stated.
- For 645 drugs (about 3–4 requests each, less with batched synonyms) this is a few minutes at ≤ 2.5 req/s; cache everything.

## 4. Are Decagon's IDs flat (`CIDm`) or stereo (`CIDs`)?

**STITCH convention** [verified, STITCH README http://stitch-db.org/download/README, reached via the redirect from http://stitch.embl.de/download/README]
- "CIDs / CID0... - this is a stereo-specific compound, and the suffix is the PubChem compound id."
- "CIDm / CID1... - this is a "flat" compound, i.e. with merged stereo-isomers. The suffix (without the leading "1") is the PubChem compound id."
- "the download files contain the prefixes CIDs/CIDm, while the API still returns CID0/CID1."
- The README also says STITCH merges salt forms and isomers, and that isomers can be looked at separately since STITCH 3. The STITCH 5 paper (Szklarczyk et al. 2016, PMC4702904) only says the chemical count is given "not including different stereoisomers"; it does not describe the ID format. The SIDER paper (Kuhn et al. 2016, PMC4702794): "In STITCH, stereoisomers and salt forms are usually merged into a common identifier … the SIDER download files contain identifiers both with and without stereochemistry."
- So the Biostars source from the first version was right about the convention.

**SIDER: flat and stereo use different numbers** [verified on `meddra_all_se.tsv.gz` (2.4 MB) and `drug_names.tsv` from https://sideeffects.embl.de/download/, README https://sideeffects.embl.de/media/download/README: "1 & 2: STITCH compound ids (flat/stereo …)"]
- 1 556 distinct (flat, stereo) pairs, 1 430 flat IDs (all `CID1…`), 1 556 stereo IDs (all `CID0…`). In 966 pairs the number is the same; in **590 pairs it differs**. 119 flat IDs have more than one stereo ID.
- Examples: ampicillin `CID100002173` / `CID000006249`; citalopram `CID100002771` / `CID000002771` and `CID000146570`; omeprazole `CID100004594` / `CID000004594` and `CID009579578`; warfarin `CID105001396` / `CID005001396`; amlodipine has two flat IDs, `CID100002162` and `CID100060496` (60496 = amlodipine besylate).
- Flat numbers can be salt or outdated records: SIDER's warfarin flat number 5001396 is warfarin **sodium** in PubChem today (current warfarin is 54678486).

**Which variant Decagon uses** [partly corrected, rest blocked]
- The Decagon paper (Zitnik et al. 2018, PMC6022705, full text via Europe PMC) names STITCH, SIDER, OFFSIDES and TWOSIDES as sources and says only that the final network was built "after linking entity vocabularies used by different databases"; it says nothing on flat vs. stereo. The SNAP page (http://snap.stanford.edu/decagon/) lists the files without an ID description. [verified: neither source states it]
- **[corrected]** The first version read `CID000002173` as "stereo / single-compound form". PubChem CID 2173 is the ampicillin anion with **no defined stereocentres** (InChIKey block 2 `UHFFFAOYSA`), and 2173 is exactly the number SIDER uses for ampicillin's **flat** ID (`CID100002173`); ampicillin's stereo ID is `CID000006249`. If `CID000002173` in Decagon is ampicillin, Decagon carries the flat PubChem number behind a `0` flag. [to verify: that `CID000002173` occurs in Decagon and is ampicillin needs the Decagon file, see *Blockers*; the ID came from the issue text, not from the file.]
- Whether all 645 IDs start with `CID0` [blocked]: `bio-decagon-combo.tar.gz` is 35 688 667 bytes (HEAD, 2026-10-28). Check once available: `cut -d, -f1,2 bio-decagon-combo.csv | tr , '\n' | grep -v STITCH | sort -u | cut -c1-4 | sort | uniq -c`.
- Practical consequence: it no longer matters much for matching, because the reverse hop (section 3) climbs from whatever CID Decagon has to an RxNorm `IN` through parent and same-connectivity steps. The flat / stereo question now matters only for reporting (which compound a Decagon node "is").
- The earlier line "Decagon strips nothing itself; its drug IDs come from TWOSIDES" is still [to verify]: the TWOSIDES paper (Tatonetti et al. 2012, PMC3382018) has no open full text in Europe PMC.

**What it means for matching**
- Enantiomer and racemate pairs (omeprazole / esomeprazole, citalopram / escitalopram) collapse at the same-connectivity step. Only use step 3 of the fallback chain when steps 1–2 found nothing, and log the matched CID and UNII so a hand check can see which form was taken.
- Salt records (warfarin sodium, amlodipine besylate) are handled by the parent step and by `PIN → IN`.

## Recommended route (summary)

**Decagon side (once, 645 drugs)**
1. STITCH ID → PubChem CID: strip `CID` and the flag digit, drop leading zeros.
2. CID → UNII (batched `compound/cid/synonyms`, keep `UNII-…`), then RxNav `rxcui.json?idtype=UNII_CODE&id=…`.
3. If nothing: parent CID (`cids_type=parent`) → step 2. If still nothing: same-connectivity CIDs of the parent → step 2. If still nothing: name match (`rxcui.json?name=…&search=2`, then `approximateTerm`).
4. Every RxCUI → `related.json?tty=IN` (and `historystatus.json` if not active) → `IN` RxCUI(s).
5. Save: STITCH ID, CID, matched CID, UNII, RxCUI, name, match step. Unmatched ones get a reason and a hand check.

**NHANES side (per distinct `RXDDRGID`)**
1. Drop `a…` / `c…` codes into the "unmappable by code" log.
2. `GN + RXDDRGID` → `findRxcuiById` (`MMSL_CODE`). On failure: `allsrc=1`, then name match on `RXDDRUG`.
3. `MIN` / `PIN` → `related.json?tty=IN`; non-active RxCUI → `historystatus.json` → `remappedConcept`.

**Join** on the `IN` RxCUI. Cache every response and record the RxNorm version (`/REST/version.json`, `05-Oct-2026` on 2026-10-28).

**Rate limits**
- RxNav: "no more than 20 requests per second per IP address", caching for 12–24 h recommended; 100 req/s with a UTS API key; RxNav-in-a-Box (Docker) for bulk use. [verified, https://lhncbc.nlm.nih.gov/RxNav/TermsofService.html]
- PubChem: 5 req/s, HTTP 503 when exceeded, watch `X-Throttling-Control`. [verified, section 3]

**Known failure cases (to log and count)**
- NHANES `a…` (not in Multum) and `c…` (unspecified, class only) codes.
- Retired Multum codes not present in the current RxNorm MMSL source.
- Decagon CIDs with no UNII on the exact, parent or same-connectivity records, and with an IUPAC-only PubChem title (name fallback fails too).
- Same-connectivity step reaching more than one `IN` (enantiomers).
- UNII that exists in PubChem but has no RxNorm concept (live: `SFV84J2559`, L-(+)-ampicillin).
- Biologics (insulins, monoclonal antibodies): no small-molecule CID; Decagon (STITCH) has none of these anyway.

## Blockers

| What | Exact request | Error / size |
|---|---|---|
| Decagon combo file (all 645 IDs; `CID0` check; is `CID000002173` ampicillin) | `curl -sSIL http://snap.stanford.edu/decagon/bio-decagon-combo.tar.gz` | `Content-Length: 35688667` (35.7 MB), over the 5 MB limit, not downloaded. Smaller official files with STITCH IDs exist on the same page but were not in the task, so not downloaded: `bio-decagon-targets.tar.gz` 61 835 B, `bio-decagon-mono.tar.gz` 1 928 066 B, `bio-decagon-targets-all.tar.gz` 441 556 B. They cover only drugs with targets / mono side effects, not necessarily all 645. |
| Share of NHANES codes that fail to map | `curl -sSIL https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/RXQ_RX_J.xpt` | `content-length: 9627680` (9.6 MB), over the limit; `RXQ_RX_I.xpt` (2015–16) is 9 629 600 B. The drug dictionary `1988/DataFiles/RXQ_DRUG.xpt` is 3 372 240 B (under the limit, not downloaded: it gives the share of distinct codes, not of prescriptions). Parsing `.xpt` also needs a SAS XPORT reader; `pandas` is not installed. |
| STITCH website over HTTPS | `curl https://stitch.embl.de/` | `curl: (60) SSL: no alternative certificate subject name matches target host name 'stitch.embl.de'`. Worked around: `http://stitch.embl.de/…` answers `301` to `http://stitch-db.org/…`, where the README was read. |
| PubChem HTML docs | `curl https://pubchem.ncbi.nlm.nih.gov/docs/programmatic-access` (also `/docs/pug-rest`, `/docs/dynamic-request-throttling`) | HTTP 200 but a Google reCAPTCHA challenge page, not the docs. Not bypassed. Worked around: the pages' own Markdown versions under `/pcfe/docs/markdown/*.md` (linked as `rel="alternate" type="text/markdown"` from the web.archive.org copy) returned the docs. |
| TWOSIDES paper full text (where Decagon's IDs come from) | `curl https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3382018/fullTextXML` | `HTTP 500` (Europe PMC lists it as not open access). |
| Decagon paper PDF | `https://arxiv.org/pdf/1802.00543` fetched (4.4 MB) | No PDF text tool available (`pdftotext`, `pdftoppm`, `pypdf` missing). Worked around: full text from Europe PMC (`PMC6022705/fullTextXML`). |

## Follow-up 2026-10-28: checks on the downloaded Decagon file

`bio-decagon-combo.tar.gz` (35,688,667 bytes, `http://snap.stanford.edu/decagon/bio-decagon-combo.tar.gz`) was downloaded to `data/raw/decagon/` in the main checkout (git-ignored) and read with a stdlib script.

- **All 645 Decagon IDs start with `CID0`.** [verified: `bio-decagon-combo.csv`, 645 distinct IDs in `STITCH 1` / `STITCH 2`]
- **The numbers behind the `CID0` prefix are mostly flat (stereo-unspecific) compound numbers, not stereo ones.** Compared with SIDER 4.1 `meddra_all_se.tsv.gz` (column 1 = flat `CID1…`, column 2 = stereo `CID0…`): 572 of the 645 numbers occur as a SIDER flat number and 380 as a stereo number. Among numbers that occur in only one of the two sets, **202 are flat-only and 10 are stereo-only**. [verified: both files, stdlib join on the numeric part] So the `CID000002173` → flat number pattern seen for ampicillin is the rule, not an exception. The README's "`CID0` = stereo" does not hold for Decagon's IDs.
- **Consequence for the reverse route:** stripping `CID` and leading zeros gives a flat PubChem CID. Flat CIDs often have no UNII (as for 2173), so the fallback chain (parent CID → same-connectivity CIDs → RxNav name match) is the normal path for many drugs, not a rare one. The 10 stereo-only numbers need a hand check.
- **Still blocked:** the share of NHANES codes that fail to map. `RXQ_RX_J.xpt` (9,627,680 bytes) is now downloaded to `data/raw/nhanes/`, but reading it needs pandas (not yet installed).

## Follow-up 2026-10-28: share of NHANES 2017–18 codes that reach an RxNorm ingredient

`RXQ_RX_J.xpt` (NHANES 2017–18) read with pandas; each distinct `RXDDRGID` looked up live with `rxcui.json?idtype=MMSL_CODE&id=GN<code>`, then `rxcui/<rxcui>/related.json?tty=IN`. Every response is cached in `data/cache/rxnav_nhanes_2017.json` (main checkout, git-ignored). [verified: live RxNav, 2026-10-28]

| Outcome | Codes | Share | Mentions | Share |
|---|---|---|---|---|
| reached IN | 632 | 94.9% | 14,214 | 99.4% |
| no RxCUI | 20 | 3.0% | 66 | 0.5% |
| `c` code (class only) | 13 | 2.0% | 19 | 0.1% |
| `a` code | 1 | 0.2% | 1 | 0.0% |

14,300 rows carry a drug code; 666 distinct codes. The most frequent failures are **pre-mixed insulin combinations** (`d04374` insulin isophane / regular, `d04510`, `d04839`), topical fluoride, penicillin (`d00116`), pancrelipase, and **`h` codes** such as `h00024` "antihypertensive agents – unspecified" (a third unmappable prefix besides `a` and `c`). Mapping failure on the NHANES side is therefore small (0.6% of mentions); coverage will be decided by the join to Decagon's 645, not by this hop. This counts only the NHANES → IN hop, not the join.
