# Preknowledge: Track B (Burden) of the polypharmacy side-effect project

**Status:** In progress — 1/5 done

## How to work through this file

<!-- Instructions for the agent. Read this section before doing anything else. -->

This is a learning map. The user is learning everything they need before starting the idea below. A topic is checked off only once the user has shown they understand it. Explain in the language the user writes in.

**Pick a topic.** Take the first open topic (`- [ ]`) unless the user names another. Re-read Destination, Why, and Already known first.

**Work it.**

1. Explain the topic briefly, at its Depth, tied to the Destination. Build on Already known and on the takeaways of finished topics.
2. Link one primary source (official docs, a standard reference, a recognised expert) to read or watch.
3. Answer questions until the user is ready.
4. Ask the **Done when** question. The user answers in their own words, from memory.
5. If the answer is shaky, re-explain the weak part and ask a different check question. Repeat until it holds up.

**Close it** only when the answer holds up **and** the user agrees it is done:

- Change `- [ ]` to `- [x]` and add a `Takeaway:` line of 1–3 sentences, phrased from the user's own answer.
- Update the Status count (skipped topics are not counted).
- Ask whether to continue with the next topic or stop.

**Keep the map honest.**

- A new topic surfaces → add it in prerequisite order and tell the user.
- A topic turns out unnecessary → mark it `- [~] skipped: <reason>`, keeping it in the list.
- Something in Not yet clear becomes sharp → turn it into a topic and remove it from Not yet clear.

**Ready check.** When every topic is `[x]` or `[~]` and Not yet clear is empty, ask 3–5 check questions mixed across the topics, answered from memory. All hold up → set Status to `Ready` and tell the user they have what they need to start. Any miss → reopen those topics (back to `- [ ]`, takeaway removed) and work them again.

## Destination

The user understands what Track B does (NHANES co-prescriptions → survey-weighted pair prevalence × the model's calibrated `p_r(a, b)` → burden ranking, novelty check, pair-level split by age and drug count; the model is a black box) well enough to make and defend their own choices when the Burden tickets are grilled: NHANES cycles (#15), the model–burden interface (#17), the burden formula (#18), and what counts as undocumented (#19), plus the map's "Not yet specified" items (RQ4 bands, failed mappings, ranking validation). Context: `PROJECT_CONTEXT.md` (`1`, `4`, `7`) and the wayfinder map, GitHub issue #1.

## Why

A personal project. The user wants to understand what they are doing so they can steer it in the right direction. Concepts and tools matter here, not medical background: keep explanations at concept level, grounded in the project's own data (NHANES drug records, Decagon's 645 drugs and 964 side effects).

## Already known

- Python data wrangling: some pandas. Has not yet shown how to write a cached, resumable, rate-limited API lookup loop (covered in topic 2).
- Burden intuition: knows prevalence × `p` is a risk weighted by how many people take the pair. Does not yet know what `p` must satisfy for that product to mean something (covered in topic 3).
- From the project context: hands-on PyTorch / PyG and GNN basics (Track A side; not needed for Track B).

## Topics

- [x] **Survey weights and pooling NHANES cycles**
  - Why: a pair's prevalence is a weighted sum over participants. This sets up the cycle choice (#15) and the age / drug-count bands (pair-level RQ4). Covers: what a weight means, why rows can't simply be counted, how weights change when cycles are pooled (incl. the 2017–Mar 2020 `P_` file), what strata and PSUs are for, subgroup analysis without subsetting, and the tools (R `survey` / Python `samplics`) at awareness level. Facts: issue #11.
  - Depth: working
  - Done when: A participant has interview weight 30,000 in 2015–16 and you pool 2015–16 with the 2017–Mar 2020 file. What does the weight mean, what happens to it in the pool, and why do you need strata and PSUs to put an error bar on a pair's prevalence?
  - Takeaway: A weight is how many Americans a participant represents. When pooling, each file's weights alone add up to the whole population, so multiply each weight by (cycle length / total years pooled), or the estimate is counted once per file. Strata and PSUs tell the software who was sampled together (similar doctors and prescribing), so the error bars aren't too narrow.

- [ ] **From NHANES drug records to Decagon IDs**
  - Why: a drug that doesn't map gets no score. This is the coverage pilot (#16), the failed-mappings decision, and how pairs outside the 645 are represented (#17). Covers: the reverse route chosen by the user (Decagon's 645: STITCH ID → PubChem CID → RxNorm ingredient, mapped once; NHANES: Multum code (`RXDDRGID`) → RxNorm, PIN → IN, MIN split into INs; join on the ingredient), where it breaks and the two failure reasons (*not in the model* vs. *mapping failed*), the four coverage numbers (drugs, mentions, participants with ≥ 2 drugs, weighted pairs), and a cached, resumable, rate-limited lookup loop. Facts: issue #12 and its route comment.
  - Depth: working
  - Done when: Walk the NHANES record "amlodipine; benazepril" to Decagon IDs with the reverse route. Name each hop, one place it can fail, which of the two failure reasons that is, and which coverage number it lowers. How do you write the lookup so a crash halfway costs nothing?

- [ ] **What the model's probability means, and the burden formula**
  - Why: this is the burden formula (#18) and what the model has to hand over (#17). Covers: `p` as a calibrated probability, prevalence × `p` as an expected count, why a good ranking isn't enough, what summing over 964 overlapping side effects counts, why every side effect counts equally (severity weights are out of scope, deferred to #20), label vs. prediction for pairs Decagon already labels, and storing each participant's full drug set so drug sets can be added later.
  - Depth: working
  - Done when: With `p` = 0.10 for nausea and 2M people taking the pair, what is 200,000, and what must be true about 0.10? Summed over all 964 side effects, what does the total count, and why could a pair rank top just from many mild side effects?

- [ ] **"Undocumented" as a matching problem**
  - Why: the novelty claim in RQ3 (#19). Covers: Decagon's own labels vs. an external reference, references that say only "these interact" (DDInter) vs. ones that name the side effect (CRESCENDDI), why TWOSIDES isn't independent evidence, and why "not found" can mean the reference simply doesn't include that drug. Facts: issue #13.
  - Depth: working
  - Done when: Your top pair has a high bleeding burden. DDInter lists it as "Major" with no side effect named, and it is not in Decagon's labels. Is it documented? Give two definitions and what each lets you claim.

- [ ] **Uncertainty and ranking stability**
  - Why: uncertainty reporting (#18) and validating the ranking (`7.1`: MEPS robustness check, positive controls). Covers: the two sources of uncertainty (survey sampling of prevalence, model `p`), bootstrap at concept level, top-k overlap between two rankings, and positive controls as a sanity check.
  - Depth: working
  - Done when: Pair A ranks #3 and pair B ranks #5. When can you say A really carries more burden than B? How do you check that your top 20 holds up when NHANES is swapped for MEPS?

## Not yet clear

None so far.
