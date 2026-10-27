# Sources for documented drug-pair interactions without a DrugBank licence

Ticket: #13 (map #1). Serves RQ3 (`1.4`) and the novelty check (`7.7`): for a top-burden **drug pair**, is the interaction "already documented"?

Researched 2026-10-27. Web access was intermittent during this session. Each claim is marked:

- **[V]** verified this session against the cited primary page or paper;
- **[K]** from prior knowledge of the source, **not re-checked this session**. Re-check before relying on it.

## Answer in one paragraph

Without a DrugBank licence, the strongest free reference is **DDInter** (curated, pair-level, severity level + mechanism text, CSV download, academic/non-commercial terms). For **side-effect-level** "documented" evidence, **CRESCENDDI** (drug pair + specific adverse event, RxNorm + MedDRA IDs) is the only clinically curated open set found, but it is small (454 drugs). **TWOSIDES / nSIDES** gives side-effect-level pairs but is the *same FAERS-mined source Decagon was built from* (`3.4.1`), so it cannot serve as independent "documented" evidence. The **ONC high-priority list** is tiny (about 15 class-level pairs) and useful only as a sanity check (`7.1`). The **NLM RxNav interaction API is gone** (discontinued around 2 January 2024) with no NLM replacement. **SIDER** has no pair data. **Drugs.com** and **KEGG** bulk use are not permitted under their terms.

## Comparison table

| Source | Licence / access | Coverage | Format | Drug IDs | Granularity | Last update |
|---|---|---|---|---|---|---|
| DrugBank (full XML) | Academic licence, CC BY-NC 4.0; downloads currently paused [V] | ~all approved drugs [K] | XML | DrugBank IDs | Pair + one-line description text [K] | 5.1.22, 2026-06-27 [V] |
| DrugBank Open Data | CC0, no licence [V] | Vocabulary + structures only, **no DDIs** [V] | CSV, SDF | DrugBank IDs, synonyms; some external IDs [K] | none | 5.1.22 [V] |
| DDInter 1.0 | "Open-access", footer "All Rights Reserved" [V]; CC BY-NC-SA 4.0 per a third-party plugin repo, not seen on DDInter's own page [V] | 1,833 drugs, 236,834 pairs [V]; downloads cover only 8 ATC groups (A, B, D, H, L, P, R, V) [V] | 8 CSVs by ATC letter [V] | DDInter IDs + names [K]; web pages link DrugBank / PubChem / ChEMBL [K] | Pair + severity (Major / Moderate / Minor / Unknown) in CSV [K]; mechanism + management text on web pages [V] | 2020-09-04 [V] |
| DDInter 2.0 | Free web access [V]; download terms not checked | 2,310 drugs, 302,516 DDI records, 8,398 mechanism/management descriptions; also drug–food, drug–disease, duplications [V] | Web; download not checked | as 1.0 [K] | as 1.0 | data 2024-05-14 (Database Commons); paper NAR 53(D1):D1356, 2024 [V] |
| CRESCENDDI | Open, Scientific Data 2022 [V]; exact licence not checked (Sci Data datasets are usually CC BY) [K] | 10,286 positive + 4,544 negative controls; 454 drugs; 179 adverse events [V] | CSV on a public repository [K] | **RxNorm** drugs, **MedDRA** events [V] | **Pair + specific adverse event** [V] | 2022 [V] |
| TWOSIDES / nSIDES | Free download from Tatonetti lab (nsides.io) [K]; licence not checked | ~59k pairs, ~1.3k events in the 2012 release [K] | CSV | RxNorm + MedDRA in nSIDES release; STITCH in Decagon's copy [K] | Pair + adverse event + PRR statistics [K] | 2012 paper; nSIDES v0.1 later [K] |
| ONC high-priority DDI list | Public, in the paper (Phansalkar et al., JAMIA 2012) [K] | ~15 class-level interactions [K] | Paper table | Drug classes / names, map by hand [K] | Pair (class) + clinical consequence text [K] | 2012 [K] |
| NLM RxNav Interaction API | **Discontinued** on or about 2024-01-02; Interactions tab removed [V, via NLM change log as quoted by search] | n/a | n/a | RxNorm (was) | Pair; data was a DrugBank subset + ONC High, no severity [V, DrugBank blog] | gone |
| SIDER 4.1 | CC BY-NC-SA 4.0 [K] | ~1.4k drugs, single-drug side effects from labels [K] | TSV | STITCH CIDs, MedDRA [K] | **Single drug only, no pairs** | 2015 [K] |
| Drugs.com | Terms forbid automated harvesting; data licensed from Cerner Multum [K] | large | web | none public | pair + severity | n/a |
| KEGG DRUG DDI | Web/API free for academic use; bulk/FTP needs subscription [K] | Japanese package inserts (contraindications, precautions) [K] | REST `rest.kegg.jp/ddi/<id>` [K] | KEGG D-numbers; links to PubChem SID, DrugBank [K] | Pair + "CI" / "P" class + mechanism note [K] | continuously [K] |

## Per source

### DrugBank: what is free without a licence
- The release page lists DrugBank **5.1.22 (2026-06-27)** and states that "All Academic DrugBank dataset downloads are temporarily paused". [V] <https://go.drugbank.com/releases/latest>
- Academic files (full XML with interactions, structures, target links) are **CC BY-NC 4.0** and need a free account plus eligibility (students, professors, research associates at an academic institution, non-commercial). [V] Same page.
- **Open Data (CC0)**: *DrugBank Vocabulary* (CSV, identifiers, names, synonyms) and *DrugBank Structures* (SDF). Neither contains drug–drug interactions. [V] Same page.
- So "free academic access" still means applying for the academic licence; without it, DrugBank DDIs are unavailable. The CC0 vocabulary is still useful as a name/synonym → DrugBank-ID bridge (DDInter and others link to DrugBank IDs).
- Granularity of the licensed DDIs: one row per pair with a templated sentence (e.g. "The risk or severity of adverse effects can be increased when ..."), no MedDRA event code. [K]

### DDInter (1.0 and 2.0)
- 1.0: 1,833 approved drugs, 236,834 associations from literature and drug labels; each has mechanism description, management, alternatives; last update 2020-09-04. [V] <https://ddinter.scbdd.com/>
- 1.0 download page: 8 CSVs named `ddinter_downloads_code_{A,B,D,H,L,P,R,V}.csv`, apparently grouped by ATC first letter. **No files for C (cardiovascular), N (nervous system), J (anti-infectives), G, M, S** are listed. [V] <https://ddinter.scbdd.com/download/>. Cardiovascular and CNS drugs dominate NHANES co-prescription, so check whether those pairs still appear (as the B/A/etc. partner) or are missing.
- Licence: no licence on the download page; footer "All Rights Reserved"; describes itself as open access. [V] A third-party BioThings plugin repo states CC BY-NC-SA 4.0. [V] <https://github.com/mnarayan1/DDInter>. The `/terms/` page could not be fetched this session.
- 2.0: Tian et al., NAR 2025, 53(D1):D1356–D1362, doi:10.1093/nar/gkae726. 2,310 drugs, 302,516 DDI records, 8,398 mechanism/management descriptions, plus 857 drug–food, 8,359 drug–disease, 6,033 therapeutic duplications; curated by pharmacists. Data last updated 2024-05-14 per Database Commons. [V] <https://ddinter2.scbdd.com/>, <https://ngdc.cncb.ac.cn/databasecommons/database/id/7473>
- CSV columns (1.0): `DDInterID_A, Drug_A, DDInterID_B, Drug_B, Level`; Level ∈ {Major, Moderate, Minor, Unknown}. [K]
- Granularity: **pair-level with severity**; the clinical effect is free text (mechanism), not a coded side effect.
- ID mapping to Decagon: DDInter ID → (DrugBank ID on drug page) → PubChem CID → STITCH; or by name via RxNorm (`4.2`). [K]

### CRESCENDDI
- Kontsioti, Maskell, Dutta, Pirmohamed, "A reference set of clinically relevant adverse drug-drug interactions", Scientific Data 9:72 (2022). [V] <https://www.nature.com/articles/s41597-022-01159-y>
- 10,286 positive and 4,544 negative controls, 454 drugs, 179 adverse events, mapped to **RxNorm** and **MedDRA**; also single-drug ADRs, indications and negative drug–event pairs. [V]
- Built by extracting and aggregating clinical resources (the paper names BNF and Micromedex-type sources [K]); intended as positive controls for DDI signal detection in spontaneous reports. [V]
- Granularity: **drug pair + specific adverse event**, the same shape as a Decagon triple (a, b, r). MedDRA → UMLS CUI mapping is needed to compare with Decagon's side-effect IDs. [K]
- Independent of FAERS mining, so usable as "documented" evidence; small drug coverage limits recall.

### TWOSIDES / nSIDES
- TWOSIDES (Tatonetti et al. 2012, `12.2`) is mined from FAERS; Decagon's combo labels are a filtered TWOSIDES (`3.4.1`). [K]
- The Tatonetti lab re-hosts OFFSIDES / TWOSIDES (nSIDES project, nsides.io) with RxNorm drug and MedDRA event IDs and PRR statistics. [K] No newer TWOSIDES release than the nSIDES v0.1 files is known. [K]
- **Circularity:** using TWOSIDES as "documented" would mostly re-find Decagon labels; it is not independent clinical documentation.

### ONC high-priority DDI list
- Phansalkar et al., "High-priority drug-drug interactions for use in electronic health records", JAMIA 2012: expert-panel list of about 15 class-level interactions (e.g. QT-prolonging agents, warfarin + NSAIDs). A companion paper lists low-priority, non-interruptive DDIs. [K]
- Granularity: drug class pair + consequence text. Use as a **sanity check** for `7.1`, not as a coverage reference.

### NLM RxNav Drug Interaction API
- NLM's RxNav change log announced that the Drug-Drug Interaction API would be discontinued **on or about 2 January 2024** and the RxNav Interactions tab removed. [V, quoted via search of <https://www.lhncbc.nlm.nih.gov/RxNav/news/news.html>; the page itself returned 403 to direct fetch]
- Its data had been a limited DrugBank subset (no severity, no management) plus ONC High. [V] <https://blog.drugbank.com/nih-discontinues-their-drug-interaction-api/> (DrugBank is an interested party.)
- No NLM replacement was found. [V for "none found"; secondary sources say none exists]
- RxNorm/RxNav **name and ID APIs** remain and are still the mapping backbone (`4.2`). [K]

### SIDER
- SIDER 4.1 (2015): single-drug side effects and indications from package inserts, STITCH CIDs + MedDRA, CC BY-NC-SA 4.0. [K]
- **No drug-pair data.** Useful only to flag that a "pair" side effect is already a known single-drug effect (`3.4.3`).

### Drugs.com
- Terms of use prohibit automated scraping / data harvesting; interaction content is licensed from Cerner Multum. [K] Manual look-up of a handful of top pairs is possible; programmatic use is not.

### KEGG DRUG
- KEGG DDI data come from Japanese package inserts, labelled contraindicated (CI) or precaution (P), with a mechanism note. Queryable per drug via `https://rest.kegg.jp/ddi/<KEGG D-number>`. [K]
- KEGG REST is free for academic use; bulk download (FTP) requires a paid subscription. [K]
- IDs: KEGG D-numbers, cross-linked to PubChem SID, ChEBI, DrugBank. [K]
- Granularity: pair-level class (CI / P), no coded side effect.

## Implications for `7.7` (facts, not decisions)

- No free source gives **side-effect-level** documentation at Decagon's scale except TWOSIDES, which is circular.
- Pair-level "documented: yes/no (+ severity)" is available from DDInter at about 1.8k–2.3k drugs, which likely covers most of Decagon's 645 drugs after mapping (not measured).
- Side-effect-level documentation is available only for CRESCENDDI's 454 drugs / 179 events.
- All licensed-DDI routes (DrugBank, DDInter, SIDER, KEGG) are non-commercial academic use.

## Not verified this session (check before use)

DDInter `/terms/` page and 2.0 download files; CRESCENDDI repository, licence and file columns; nSIDES file layout and licence; ONC list contents; SIDER licence; Drugs.com terms wording; KEGG API terms. Fetch attempts to these failed for tool reasons, not because the pages are missing.

## Sources

- DrugBank releases: <https://go.drugbank.com/releases/latest>
- DrugBank blog on RxNav API: <https://blog.drugbank.com/nih-discontinues-their-drug-interaction-api/>
- NLM RxNav news: <https://www.lhncbc.nlm.nih.gov/RxNav/news/news.html>
- DDInter 1.0: <https://ddinter.scbdd.com/>, <https://ddinter.scbdd.com/download/>; Xiong et al., NAR 2022, 50(D1):D1200, doi:10.1093/nar/gkab880
- DDInter 2.0: <https://ddinter2.scbdd.com/>; Tian et al., NAR 2025, 53(D1):D1356, doi:10.1093/nar/gkae726; <https://ngdc.cncb.ac.cn/databasecommons/database/id/7473>
- DDInter BioThings plugin (licence claim): <https://github.com/mnarayan1/DDInter>
- CRESCENDDI: Kontsioti et al., Sci Data 2022, 9:72, <https://www.nature.com/articles/s41597-022-01159-y>
- Phansalkar et al., JAMIA 2012 (ONC high-priority DDIs)
- Tatonetti et al., Sci Transl Med 2012 (TWOSIDES); nSIDES: <https://nsides.io>
- SIDER: <http://sideeffects.embl.de>
- KEGG REST: <https://www.kegg.jp/kegg/rest/keggapi.html>
