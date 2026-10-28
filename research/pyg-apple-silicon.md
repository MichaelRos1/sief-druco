# PyG on Apple Silicon (MPS): feasibility of the Decagon encoders

Research for issue #3 (map #1). Researched 2026-10-27; second pass 2026-10-28 (items from the
PARTIAL comment: `scatter_reduce` amax on MPS, memory/time figures, dense RGCN rewrite).
Point IDs (`3.2`, `5.3`, ...) refer to `PROJECT_CONTEXT.md`.

**Machine checked:** Apple M4, 24 GB unified memory, macOS 15.5 (`sysctl`, `sw_vers`).
**Software state at time of writing:** PyTorch 2.14.1 is the latest GitHub release (2026-09-30) [S14];
PyPI already lists `torch` 2.15.0 (uploaded 2026-10-28) but its git tag could not be read (see
Blockers). PyG 2.8.1 is the latest release [S2].

**Status of this file: still PARTIAL.** The source-level questions are now settled from primary
sources. Nothing here has been run: this Mac has no `torch` / `torch_geometric` installed and the
Decagon data is not downloaded (section 8).

Legend:
- **[Sn]** = cited primary source (list at the end). Source links are pinned to commits/tags:
  PyG master at `ea2c5f6` (2026-10-21), PyTorch tags as named.
- **[verified: Sn]** = read in the cited source.
- **[derived]** (also written **(est.)** in older text) = my own arithmetic from the PyG source
  (which tensors are allocated) and the Decagon counts in `3.1`/`3.2`. **Not a measurement.**
- **[unverified]** = could not be confirmed from a primary source.

---

## 1. Short answer

| Encoder (`5.3`) | Runs on MPS? | Memory, 2 layers, d = 64 [derived] | Locally feasible on the M4 / 24 GB? |
|---|---|---|---|
| Tensor factorisation (no encoder) | Yes. Only `nn.Embedding`, elementwise ops and matmul. | < 0.5 GB per 64k batch | **Yes**, easily. Seconds per epoch (est.). |
| `to_hetero` GraphSAGE, collapsed drug–drug edges | Yes: sum/mean aggregation uses plain `scatter_add_` / `index_select` [verified: S3]. | Encoder **< 1 GB** (was 1–3 GB) | **Yes.** |
| `to_hetero` GAT, collapsed drug–drug edges | Yes on PyTorch ≥ 2.6: the softmax's scatter max is `scatter_reduce_(reduce='amax')`, which has an MPS kernel since 2.6.0 and a hand-written Metal kernel since 2.13.0 [verified: S4, S3, S15–S18]. Keep `edge_index` a plain `Tensor` (section 2). | **≈ 1–1.5 GB** (was 3–4 GB) | **Yes** by source; not yet run. |
| HGT, collapsed drug–drug edges (4 edge types) | Same softmax path as GAT [verified: S19]; type-specific linears run as a Python loop over heads × edge types without `pyg-lib` [verified: S6, S19]. | **≈ 3–3.5 GB** | **Yes** by source, if collapsed. With all 967 edge types: ≈ 12–15 GB and a 3,868-step loop (heads = 4). Not practical. |
| `RGCNConv`, 967 relations, basis decomposition | Runs with core ops, but loops over all 967 relations per layer; `segment_matmul` path is disabled when `num_bases` is set [verified: S7]. | **≈ 10 GB** (confirmed) | **Marginal.** Slow (est.). Better with the rewrite in 5.1. |
| `FastRGCNConv`, same setup | Same ops, but gathers one weight matrix per edge, and its mean normalisation builds a one-hot [E, 967] tensor [verified: S7]. | **≈ 108 GB** per layer for the weights, plus ≈ 51 GB transient for the one-hot normalisation | **No.** |

Bottom line: tensor factorisation, GraphSAGE, GAT and HGT on the collapsed graph run on MPS
**according to the source code** (PyTorch ≥ 2.6; ≥ 2.13 for the native Metal scatter kernels). Stock
`RGCNConv` with 967 relations is borderline; `FastRGCNConv` is not feasible. The dense rewrite of the
RGCN drug–drug block (5.1) matches `RGCNConv`'s algebra, with two caveats on how it is combined with
the protein relations. All of this still needs the run in section 6.

---

## 2. How PyG executes on MPS

- **Core PyG needs no extension libraries.** "From PyG 2.3 onwards, you can install and use PyG
  without any external library required except for PyTorch" [S8]. In that mode, message passing is
  built from core PyTorch ops, so it runs on whatever device the tensors are on, including MPS.
- **PyG's maintainer said this early on.** Core ops should get GPU acceleration on M1 for layers such
  as `GCNConv` / `SAGEConv` that do not rely on `torch-scatter` / `torch-sparse`. MPS support for
  PyG's own external operators is "of lower priority" [S9, rusty1s, 2022-06-09].
- **Which core op each aggregation uses** [verified: S3 `_scatter.py` L62–138]:
  - `sum` / `add` / `mean`: `scatter_add_` (L67–80; `mean` divides by a `clamp(min=1)` count,
    L72–80), and `scatter_` for `any`. There are no device checks.
  - `min` / `max`: `src.new_zeros(size).scatter_reduce_(dim, index, src, reduce='amax'|'amin',
    include_self=False)` (L96–100). The `torch_scatter` path (L115) is taken only if
    `torch-scatter` is installed **and** `src.is_cuda` **and** `src.requires_grad` (L85–87). On MPS
    the core PyTorch kernel is always used.
  - `mul`: `scatter_reduce_(reduce='prod')` (L119–133).
  - The only device branch in the file is `src.is_cuda`; there is **no MPS-specific code path**.
- **How the attention softmax calls it** [verified: S4 `_softmax.py`]:
  - GAT calls `softmax(alpha, index, ptr, dim_size)` in `edge_update` [verified: S20 L404]; HGT
    calls `softmax(alpha, index, ptr, size_i)` in `message` [verified: S19 L230].
  - With `index` and no `ptr`, `softmax` takes `scatter(src.detach(), index, dim, dim_size=N,
    reduce='max')` (L82–88), then `exp`, then a `scatter(..., reduce='sum')`. Because the input is
    detached, `requires_grad` is False, so this is always `scatter_reduce_(amax)` even when
    `torch-scatter` is installed [verified: S3 L85–87]. No backward of the max is needed.
  - **`ptr` is `None` for a plain `Tensor` `edge_index`** [verified: S21 `message_passing.py`
    L379–386]. It is set only for a sorted `EdgeIndex` (L387–390) or a `torch_sparse.SparseTensor`
    (L392–400).
  - **Caveat (new):** with `ptr` set, `softmax` uses `segment(...)` instead (L60–81). Without
    `torch-scatter` this calls `torch._segment_reduce` [verified: S22 L27–28, L37–50], whose
    dispatch in PyTorch 2.14 is `CPU, CUDA, XPU` only, no MPS [verified: S15b L14868–14872]; PyTorch's
    own MPS test list marks `_segment_reduce` as unimplemented [verified: S18 L368, L586–587]. With
    `torch-scatter` it calls `torch_scatter.segment_csr` (L34), a CPU/CUDA kernel. **So on MPS, pass
    `edge_index` as a plain `torch.Tensor`, not a sorted `EdgeIndex` or `SparseTensor`.**
    Ordinary aggregations are less exposed: `Aggregation.reduce` takes `segment` only when `index`
    is missing or the aggregation is deterministic, otherwise `scatter` [verified: S24
    `aggr/base.py` L172–185].
- **PyTorch MPS kernel for `scatter_reduce_(amax)` (settled):**
  - `scatter_reduce_` with `include_self` is the ATen op `scatter_reduce.two`, implemented by the
    structured kernel `scatter_reduce_two` [verified: S15].
  - **PyTorch ≤ 2.5: no MPS kernel.** In v2.5.0 the dispatch is `CPU, CUDA: scatter_reduce_two`
    [verified: S15a L8250]. On MPS this raised `NotImplementedError` or, with
    `PYTORCH_ENABLE_MPS_FALLBACK=1`, ran on the CPU.
  - **PyTorch 2.6.0: first MPS kernel.** PR #141948 "[MPS] Add scatter_reduce.two" (malfet, commit
    `d648133`, 2024-12-04) added `MPS` to the dispatch and a `scatter_reduce_two_mps_kernel` that
    redispatches to the existing MPSGraph `scatterAlongAxis` code with mode `amax` [verified: S16,
    S15c L8291, S17a L343–367]. The commit is contained in tag v2.6.0 and not in v2.5.0 (GitHub
    compare API). Limits named in the PR: min/max not supported for int64; int32 only on macOS 15+.
    float32 (our case) is supported.
  - **PyTorch 2.13.0: hand-written Metal kernel.** PR #184028 "[MPS] Migrate scatter/gather ops to
    Metal with async bounds checking" (commit `9b7a45e`, 2026-05-18, in tag v2.13.0, not in v2.12.0)
    moved the scatter family "including reduce parts" to native Metal; amin/amax use an atomic CAS
    loop [verified: S17]. In v2.14.0, `scatter_reduce_two_stub` is registered to
    `scatter_reduce_two_mps_kernel` → `scatter_reduce_dispatch` → `scatter_reduce_metal`
    [verified: S17b L189, L439–473, L488–501], and the Metal library instantiates amax for `float`
    [verified: S17c L486–495, L519]. The 2.13 blog's "scatter/gather with bounds checking" [S10]
    refers to this PR.
  - `include_self=False` is handled before the kernel, device-independently: the shared
    `scatter_impl` fills the touched rows with `-inf` for max (`scatter_reduce_exclude_self_helper`)
    [verified: S15d L2232–2263, L2133–2160]. The CPU-only "expanded index" fast path returns early
    for non-CPU tensors [verified: S15d L2037], so MPS always goes through the stub.
  - **Preliminary claim corrected:** the first pass said the 2.13 notes confirm scatter/gather "but
    not the reduce variants" and that a CPU fallback might be needed. In fact `amax` has been native
    on MPS since 2.6.0 (MPSGraph) and is a native Metal kernel since 2.13.0. Use PyTorch ≥ 2.13.
  - What the source cannot tell: speed, and whether the float atomics behave as expected on this
    M4. The smoke test in section 6 checks that it runs.
- **CPU fallback for missing ops.** `PYTORCH_ENABLE_MPS_FALLBACK=1` makes PyTorch "fall back to the
  CPU for operations that MPS doesn't support" [S5]. Without it, a missing op raises
  `NotImplementedError ... is not currently implemented for the MPS device`. An example is
  `aten::_unique2` in PyG's TGN example [S11].
- **Memory limits on MPS** [S5]:
  - `PYTORCH_MPS_HIGH_WATERMARK_RATIO` (default 1.7) is a hard cap, as a multiple of the device's
    recommended maximum working set.
  - `PYTORCH_MPS_LOW_WATERMARK_RATIO` (default 1.4 on unified memory) is a soft limit.
  - On a 24 GB machine the working set is shared with the CPU side; query it with
    `torch.mps.recommended_max_memory()` [verified: S23 L130].
- **Open PyG reports.** One user reports steadily growing memory and slowdowns with a GAT model on
  MPS (PyG 9754, open, no maintainer reply) [S12]. There is no official MPS support statement in PyG
  2.8.x release notes [S2].

## 3. Extension libraries on macOS arm64

| Library | macOS wheels (PyTorch 2.8 index) [S13] | Devices | Notes |
|---|---|---|---|
| `torch-scatter` 2.1.2 | `universal2` (includes arm64) | CPU and CUDA kernels only; no Metal kernels (unverified for the latest release) | PyG ignores it on non-CUDA tensors for min/max/mul [S3]. Of little use on MPS. |
| `torch-sparse` 0.6.18 | `universal2` | CPU and CUDA (unverified) | Only needed for `SparseTensor`. PyG 2.6 announced that `EdgeIndex`/`Index` "will fully deprecate the usage of `SparseTensor`" [S2]. |
| `pyg-lib` 0.6.0 | native `macosx_11_0_arm64` (0.4/0.5: `universal2`, macOS 14) | CPU and CUDA | Provides `segment_matmul` / `grouped_matmul` (heterogeneous ops) and sampling [S8]. From PyG 2.8.0 it also replaces `torch-cluster` / `torch-spline-conv` (requires `pyg-lib==0.7.0`) [S2]. |
| `torch-cluster`, `torch-spline-conv` | `universal2` (old versions) | — | Deprecated in PyG 2.8.0 [S2]. None of the encoders needs them. |

- **Install note.** Wheels are published per PyTorch version on `data.pyg.org` [S8]. Conda packages
  stop at PyTorch 2.5 [S8].
- **Consequence.** On MPS none of these libraries accelerates anything:
  - Their kernels run on CPU tensors only (unverified for each latest release).
  - PyG's fast paths that use them (`segment_matmul` in `RGCNConv` and `HeteroLinear`) are skipped
    or would fail for MPS tensors.
  - The minimal install (PyTorch + `torch_geometric` only) is the right baseline for MPS. Optional
    `pyg-lib` helps only if message passing runs on the **CPU**.

## 4. Per-encoder analysis (Decagon scale)

Every number in this section is **[derived]**: it counts the tensors that the cited PyG code
allocates, at the sizes below. None is measured. Times additionally depend on the hardware
assumptions in 4.1, which are themselves unverified.

### 4.1 Assumptions behind the estimates

**Graph sizes:**

| Quantity | Value |
|---|---|
| Nodes N | 645 drugs + 19,081 proteins = **19,726** |
| PPI edges, both directions | 1,431,224 |
| Target edges, + reverse | 37,380 |
| Combo triples after the 964 filter | ≈ 4.6M |

**Split:**
- 80 % of pairs go to training; within training, 70 % are message edges and 30 % are supervision
  edges (`RandomLinkSplit(disjoint_train_ratio=0.3)`, `5.1`).
- Message graph: ≈ 2.58M typed triples, which is **5.15M directed typed drug–drug edges**. These
  collapse to **≈ 71k directed drug pairs**.
- Message edges in total: E ≈ 1.43M + 0.037M + 5.15M = **6.62M** typed, or
  E_c ≈ 1.43M + 0.037M + 0.071M = **1.54M** collapsed.
- Supervision per epoch: ≈ 1.10M positives + 1.10M degree-matched negatives = **2.2M scored
  triples**.
- Tensor factorisation has no message graph, so all training triples are supervision:
  ≈ 7.4M per epoch.

(Second pass: the split arithmetic above was re-checked and is consistent.)

**Model:**
- d = 64, float32, 2 encoder layers; 4 attention heads for GAT/HGT.
- Decagon decoder z_iᵀ D_r R D_r z_j costs ≈ 2d² + 4d ≈ 8.4 kFLOP per triple forward, ≈ 25 kFLOP
  forward + backward. It saves about six d-vectors per triple for backward (z_i, z_j, gathered D_r
  and three intermediates): ≈ 1.5 KB per triple.

**Hardware (unverified):**
- M4 memory bandwidth 120 GB/s. Second pass could not re-confirm this: Apple's current spec pages
  list M5 models, and the M4 newsroom page gives no GB/s figure (see Blockers).
- FP32 GPU throughput: Apple does not publish it; ~4 TFLOPS peak is a commonly quoted figure.
- I assume 1 TFLOPS effective in the time estimates.

### 4.2 Tensor factorisation (no encoder)

- **Ops:** `Embedding`, `index_select`, elementwise, `matmul`, `BCEWithLogitsLoss`. All are standard
  MPS ops.
- **Parameters:** 645·64 + 964·64 + 64² ≈ 0.11M.
- **Time [derived]:** ≈ 7.4M × 25 kFLOP ≈ 0.18 TFLOP per epoch, so well under a second of compute
  at the assumed throughput. With ~64k-triple batches that is about 115 steps; per-step overhead
  will dominate. "A few seconds per epoch" is a guess, not a measurement.
- **Memory [derived]:** ≈ 1.5 KB of activations per triple, so ≈ 100 MB per 64k batch.
- **Verdict:** Feasible on MPS; the cheapest model.

### 4.3 GraphSAGE via `to_hetero`, drug–drug edges collapsed to one type

- **Ops:** gather (`index_select`), `scatter_add_` / mean, `Linear`. All core ops [verified: S3].
- **What is allocated per layer [derived]:**
  - Transient gathered messages x_j: E_c × d × 4 B ≈ 1.54M × 256 B ≈ **0.39 GB**. They are not kept
    for backward: the backward of `index_select` and `scatter_add_` needs only the index.
  - Kept for backward: the edge indices (≈ 25 MB) and a few [N, d] tensors (≈ 5 MB each).
- **Memory [derived]:** encoder peak **< 1 GB** for 2 layers. **Correction:** the first pass said
  1–3 GB; that overstated it, because the per-edge messages are not saved.
- **Time:** a single encoder pass moves ~1–2 GB of memory, ~0.1–0.5 s with overhead (est.).
- **Per epoch, two ways to train:**
  - Full-batch: one encoder pass plus the decoder over 2.2M triples in chunks, then one backward.
    Decoder activations ≈ 2.2M × 1.5 KB ≈ 3.4 GB [derived]. About 1–3 s (est.).
  - Minibatch: re-run the encoder per 64k batch, about 34 passes, so ~5–20 s (est.).
- **Verdict:** Feasible.

### 4.4 GAT via `to_hetero`, collapsed

- **Ops:** the same as GraphSAGE, plus `softmax` with a scatter max (`scatter_reduce_(amax)`), which
  is native on MPS from PyTorch 2.6 [verified: section 2].
- **What is allocated per edge and layer** (`GATConv(64, 16, heads=4)`; projection before
  propagation, `message = alpha.unsqueeze(-1) * x_j`) [verified: S20 L297–361, L387–409]:
  - Kept for backward: x_j [E, 64] (256 B) and about six [E, 4] attention tensors (~100 B), plus
    indices: ≈ **370 B per edge** [derived].
  - Transient: the message output [E, 64] (256 B).
- **Memory [derived]:** ≈ 0.57 GB kept per layer, ≈ 1.1 GB for 2 layers, peak ≈ **1–1.5 GB**.
  **Correction:** the first pass assumed ~1 KB per edge and got 3–4 GB; that overstated it.
- **Time:** about 2–3× GraphSAGE (est.).
- **Practical caveats:**
  - `GATConv` on bipartite edge types (drug→protein) needs `add_self_loops=False` under `to_hetero`.
  - Pass `edge_index` as a plain `Tensor` (section 2), so the softmax uses `scatter_reduce_`.
- **Verdict:** Feasible by source; the smoke test confirms it runs.

### 4.5 HGT, collapsed (edge types: ppi, targets, rev_targets, combo)

- **Ops** [verified: S19]:
  - K/Q/V projections per node type (`HeteroDictLinear`).
  - Relation transforms `k_rel` / `v_rel` are `HeteroLinear` layers with
    `num_types = heads × len(edge_types)` (L79–83), applied to **source nodes per edge type**, not
    per edge (`_construct_src_node_feat`, L119–157).
  - One `propagate` over a single bipartite edge index built from all edge types (L195–199).
  - Attention softmax (scatter max) as in GAT (L230).
- **Without `pyg-lib`:** `HeteroLinear` falls back to `forward_naive`, a Python loop
  `out[start:end] = x[start:end] @ self.weight[i]` over all `num_types` types [verified: S6
  L248–252]. **Correction:** the loop runs over heads × edge types, not edge types: 16 iterations
  collapsed, **3,868** with 967 edge types at heads = 4 (the first pass said ~967).
- **What is allocated per edge and layer [derived]:**
  - Kept for backward: gathered k_j, q_i, v_j, each [E, 64] (3 × 256 B), about six [E, 4]
    tensors, indices: ≈ **0.9 KB per edge**.
  - Transient: `q_i * k_j` and the message output, each [E, 64]: another ≈ 0.5 KB.
- **Memory [derived]:**
  - Collapsed (E_c ≈ 1.54M): ≈ 1.35 GB kept per layer, ≈ 2.7 GB for 2 layers, peak ≈ **3–3.5 GB**.
    Consistent with the first pass's 3–5 GB.
  - All 967 edge types (E ≈ 6.62M): ≈ 5.8 GB kept per layer, ≈ 11.6 GB for 2 layers, peak
    ≈ **12–15 GB**. Consistent with the first pass's ~13 GB. The per-type source features add only
    ≈ 0.17 GB per tensor (964 × 645 + 38.8k ≈ 0.66M rows).
- **Verdict:** Feasible only with collapsed edge types.

### 4.6 `RGCNConv`, 967 relations, `num_bases = 30`, on `data.to_homogeneous()`

How the PyG layer runs in this configuration [verified: S7]:

- **Defaults:** `aggr='mean'`, `root_weight=True`, `bias=True` (L99–100, plus `bias`).
- **Bases.** `weight = (comp @ basis.view(B, -1)).view(R, in, out)` once per forward (L205):
  [967, 64, 64] ≈ 16 MB.
- **`segment_matmul` is off.** It needs `pyg-lib`, no compilation, **`num_bases is None`**,
  float input and a plain `Tensor` `edge_index` (L242–246). With bases it is always off.
- **Loop over relations** (L257–270). For each of the 967 relations:
  - `tmp = edge_index[:, edge_type == i]` (`masked_edge_index`, L23–26);
  - `h = propagate(tmp, x=x_l, size=(N, N))`, a per-relation mean over N = 19,726 rows;
  - `out = out + h @ weight[i]`.
- Then `out + x_r @ root + bias` once (L275–280).
- **The docs warn about this pattern.** `RGCNConv` "can use dynamic shapes" and "may result in low
  GPU utilization in case the graph has a large number of relations".

Cost [derived], with E ≈ 6.62M edges in the homogeneous graph:

- **Masking traffic.**
  - Per iteration: read `edge_type` (int64, 53 MB), write and read the bool mask (2 × 6.6 MB),
    read `edge_index` (106 MB): ≈ 172 MB. Over 967 iterations: ≈ **166 GB per layer forward**,
    ≈ 1.4 s at 120 GB/s.
  - **Correction:** the first pass counted only `edge_index` (≈ 100 GB, 0.9 s).
  - Boolean-mask indexing has a data-dependent output size: one device synchronisation per
    relation, 967 per layer.
- **Kept for backward (confirmed).**
  - `h @ weight[i]` saves `h` of shape [N, d] for the weight gradient, for every relation, even
    though 964 of the relations only reach the 645 drugs (the other rows are zero).
  - 967 × 19,726 × 64 × 4 B ≈ **4.9 GB per layer, ≈ 9.8 GB for 2 layers** (d = 32: ≈ 4.9 GB).
  - Plus the masked edge indices (≈ 0.1 GB per layer) and per-relation mean counts
    (967 × N × 4 B ≈ 0.08 GB per layer).
- **Featureless variant.** If `x` is a node-index tensor instead of float features, the loop takes
  `propagate(tmp, x=weight[i, x_l])` (L259–264) and the basis-expanded weight becomes
  [967, N, d] ≈ 4.9 GB per layer instead.
- **Time (est.).** About 2–10 s per encoder forward + backward, dominated by masking traffic and
  the 967 synchronisations. Per epoch that is a few seconds to tens of seconds full-batch, or
  minutes with per-batch encoder passes.
- **Verdict:** Runs on MPS with core ops, but uses about 40–60 % of the realistic MPS budget on
  24 GB and is slow. Marginal.

### 4.7 `FastRGCNConv`, same setup

- **How it works** [verified: S7 L314–374]: one `propagate` over all edges with no loop over
  relations. Its `message` does `torch.bmm(x_j.unsqueeze(-2), weight[edge_type])` (L361). Its
  `aggregate` with `aggr='mean'` computes a per-relation degree via
  `one_hot(edge_type, num_relations)` → `scatter(...)[index]` (L368–371). The docs warn it "may
  consume a large amount of memory to compensate" [S7].
- **Memory [derived]:**
  - The gathered weight tensor `weight[edge_type]` is [E, 64, 64]: 6.62M × 64 × 64 × 4 B
    ≈ **108 GB per layer**, and `bmm` keeps it for backward (confirmed). Drug–drug edges alone:
    ≈ 84 GB.
  - **Missed in the first pass:** the mean normalisation allocates `one_hot` [E, 967] float32
    ≈ 25.6 GB and then `scatter(...)[index]`, another [E, 967] ≈ 25.6 GB, both alive at once:
    ≈ **51 GB transient, independent of d**.
  - **Correction:** the first pass said "even at d = 16: ≈ 7 GB per layer". The weight gather is
    ≈ 6.8 GB at d = 16, but the one-hot normalisation still needs ≈ 51 GB, so d does not rescue it
    (only `aggr='add'` avoids the one-hot).
- **Verdict:** Infeasible locally, and on most single GPUs.

## 5. Workarounds

### 5.1 Rewrite RGCN's drug–drug block as dense drug-only algebra (my derivation, checked against the source)

The 964 side-effect relations connect only the 645 drugs. With basis decomposition
W_r = Σ_b a_rb V_b, the relational sum for drug i is a reordering of the same sum:

Σ_r Σ_{j ∈ N_r(i)} (1/|N_r(i)|) x_j W_r = Σ_b [ (Σ_r a_rb Â_r) X_drug ]_i V_b

Here Â_r is the row-normalised 645 × 645 adjacency of side effect r, so the per-relation mean
normalisation is folded in. (The first pass wrote the left side as W_r x_j, column-vector style,
and the right side row-vector style; PyG uses row vectors, `h @ weight[i]`. The algebra is the same.)

**Implementation:**
1. Precompute `A = stack(Â_r)` as a constant dense tensor [964, 645, 645]. This is ≈ 1.6 GB in
   float32, or 0.8 GB in float16 [derived].
2. Each layer:
   - `M = einsum('rb,rij->bij', comp, A)` gives [30, 645, 645] (≈ 50 MB). That is ≈ 24 GFLOP
     forward and about as much again for the gradient with respect to `comp` [derived].
   - `out_drug = Σ_b M_b @ X_drug @ V_b` (≈ 1.6 GFLOP) [derived].
3. Handle the 3 protein relations (PPI, targets, reverse targets) separately and add the results —
   see the caveats below for how.

**Check against the `RGCNConv` source** [S7, S3]:

| Aspect | `RGCNConv` (source) | Rewrite | Match? |
|---|---|---|---|
| Basis decomposition | `weight = comp @ basis` → W_r = Σ_b comp[r, b] · basis[b], shape [in, out] (L205) | Σ_b a_rb V_b with `a = comp`, `V = weight` | **Yes** |
| Per-relation message | `h = propagate(edge_index[:, edge_type == i], x)`, then `h @ weight[i]` (L257–270) | `(Â_r X) W_r`, reordered into Σ_b M_b X V_b | **Yes** (linear reordering) |
| Normalisation | `aggr='mean'`: `MeanAggregation` → `scatter(..., reduce='mean')` (S24), i.e. sum of x_j over r-neighbours divided by the count, count `clamp(min=1)` (S3 L72–80), i.e. 1/\|N_r(i)\| per relation, zero for nodes without r-neighbours | Row-normalised Â_r, zero rows stay zero | **Yes**, if Â_r counts duplicate edges (mean counts multi-edges) |
| Direction | Default flow `source_to_target`: node i = `edge_index[1]` aggregates from j = `edge_index[0]` | Â_r[i, j] must be built with row = target, column = source | **Yes**, if built that way; a transposed Â_r is wrong unless the drug–drug graph is exactly symmetric |
| Root weight and bias | Added **once per layer**: `out + x_r @ root + bias` (L275–280) | Not in the formula | **Differs:** add `x @ root + bias` once. If the protein relations use a second `RGCNConv`, give that one `root_weight=False, bias=False` (or vice versa); a plain "add the results" counts root and bias twice |
| Shared bases across all 967 relations | One `comp` [967, 30] and one basis set for all relations, including PPI / targets / reverse targets | Step 3 with a separate `RGCNConv(num_relations=3)` gives the protein relations their own weights | **Differs:** identical only if the 3 protein W_r are also built from the same V_b (rows 964–966 of `comp`). Otherwise it is a slightly different model (more parameters), not "mathematically identical" as the first pass said |
| Input type | Float features: `h @ weight[i]`; node-index input: `weight[i, x_l]` path (L259–264) | Assumes float `X_drug` | Featureless input needs `X = I` (or an embedding) explicitly |
| Gradients | Through `comp`, `basis`, `root`, `bias` | `A` is a constant, so gradients reach `comp` via the einsum and `V` via the matmul | **Yes** |

**Properties:**
- Equivalent to `RGCNConv` with bases on the 964 drug–drug relations, provided root/bias and the
  shared bases are handled as in the table.
- Uses only dense `einsum`/`matmul`, the best-supported work on MPS.
- ≈ 1.7 GB of memory (A + M + gradients) and, at the assumed 1 TFLOPS, ≈ 50 ms per layer forward +
  backward [derived].
- Removes the 967-iteration loop and the ~10 GB of saved activations.
- **Not yet tested.** The equivalence script in section 6 compares it with `RGCNConv` on a small
  random graph (forward output and the gradient of `comp`).

### 5.2 Message passing on the CPU, decoder (and TF) on MPS

- **How:** autograd tracks `.to()` across devices. `z = encoder_cpu(...)` then `z.to('mps')`, with
  the loss on MPS and `loss.backward()`, works.
- **Transfer cost:** only z_drug [645, 64] ≈ 165 KB per step.
- **Why it can help:**
  - PyG's own comment says `torch.scatter_reduce` is faster on CPU than `torch-scatter` [S3].
  - On the CPU, `pyg-lib` CPU kernels (e.g. `segment_matmul` for `HeteroLinear`) become usable [S6].
- **Costs:**
  - The CPU and GPU share the same unified memory bandwidth.
  - Bandwidth-bound scatter work may therefore not speed up much in either place.

### 5.3 `PYTORCH_ENABLE_MPS_FALLBACK=1`

- **Effect:** any op MPS lacks (e.g. `aten::_unique2` [S11], `_segment_reduce` [S15b]) runs on the
  CPU instead of raising an error [S5]. Slower, but correct.
- **How to set it:** in the environment before starting Python.
- **Debugging tip:** leave it **off** while testing, so missing ops show up as errors.

### 5.4 Subsample side effects

- **Effect on RGCN:** loop length, saved activations and masked-edge traffic all scale linearly with
  the number of relations. Top 100 side effects: ≈ 0.5 GB per layer instead of 4.9 GB [derived].
- **Cost:** results are then no longer comparable with Decagon's 964-type numbers (`6.2`). Use it
  for development and debugging, not for the final comparison.

### 5.5 Smaller d, full-batch training, mixed precision

- **Smaller d:** d = 32 halves the per-edge and per-node activation figures above (not the
  `FastRGCNConv` one-hot).
- **Full-batch training:** one encoder pass per epoch, with the decoder loss accumulated over chunks
  of supervision triples. This removes the dominant "encoder re-run per batch" factor.
- **Mixed precision:** float16 for the constant `A` in 5.1 halves its memory. float64 should be
  avoided entirely. A float64 multi-hot drug feature matrix from pandas must be cast to float32
  before `.to('mps')`. The MPS backend does not support float64 (unverified: not checked against a
  primary source here).

## 6. Tests to run on the MacBook (not run yet)

These need PyTorch and PyG, which are not installed on this Mac (section 8). Use PyTorch ≥ 2.13 so
the Metal scatter kernels are used; 2.14.1 is the newest release whose source was checked.

**Install** (into a fresh virtual environment; `torch` 2.14.1 has a `cp313` macOS arm64 wheel on
PyPI, ≈ 127 MB; PyG's install command is `pip install torch_geometric` [S8]):

```bash
cd /Users/mr/Aspace/sief-druco && python3 -m venv .venv
source .venv/bin/activate
python -m pip install "torch==2.14.1" "torch_geometric==2.8.1"
```

### 6.1 MPS smoke test

Run once with **no** `PYTORCH_ENABLE_MPS_FALLBACK` set (`env -u PYTORCH_ENABLE_MPS_FALLBACK python
smoke_mps.py`). Any `NotImplementedError` names the missing op. Expected from the source: all `OK`
except `softmax_ptr` (`_segment_reduce` has no MPS kernel).

```python
# smoke_mps.py
import torch, torch_geometric
from torch_geometric.nn import SAGEConv, GATConv, HGTConv, RGCNConv
from torch_geometric.utils import scatter, softmax

print('torch', torch.__version__, '| pyg', torch_geometric.__version__,
      '| mps', torch.backends.mps.is_available())
dev = torch.device('mps')
N, E, d = 2000, 20000, 64
x = torch.randn(N, d, device=dev)
ei = torch.randint(0, N, (2, E), device=dev)
et = torch.randint(0, 967, (E,), device=dev)
idx = ei[1]
ptr = torch.cat([torch.zeros(1, dtype=torch.long),
                 torch.bincount(idx.cpu(), minlength=N).cumsum(0)]).to(dev)

def hgt():
    from torch_geometric.data import HeteroData
    h = HeteroData()
    h['drug'].x = torch.randn(64, d); h['protein'].x = torch.randn(500, d)
    h['drug', 'combo', 'drug'].edge_index = torch.randint(0, 64, (2, 800))
    h['protein', 'ppi', 'protein'].edge_index = torch.randint(0, 500, (2, 3000))
    h['drug', 'targets', 'protein'].edge_index = torch.stack(
        [torch.randint(0, 64, (300,)), torch.randint(0, 500, (300,))])
    h['protein', 'rev_targets', 'drug'].edge_index = h['drug', 'targets', 'protein'].edge_index.flip(0)
    h = h.to(dev)
    conv = HGTConv(d, d, h.metadata(), heads=4).to(dev)
    return conv(h.x_dict, h.edge_index_dict)['drug']

tests = {
    'scatter_reduce_amax': lambda: torch.zeros(N, d, device=dev).scatter_reduce_(
        0, idx.view(-1, 1).expand(E, d), x[ei[0]], reduce='amax', include_self=False),
    'scatter_sum': lambda: scatter(x[ei[0]], idx, 0, dim_size=N, reduce='sum'),
    'scatter_max': lambda: scatter(x[ei[0]], idx, 0, dim_size=N, reduce='max'),
    'softmax_index': lambda: softmax(torch.randn(E, device=dev), idx, num_nodes=N),
    'softmax_ptr': lambda: softmax(torch.randn(E, device=dev), None, ptr, num_nodes=N),
    'SAGEConv': lambda: SAGEConv(d, d).to(dev)(x, ei),
    'GATConv': lambda: GATConv(d, d // 4, heads=4).to(dev)(x, ei),
    'HGTConv': hgt,
    'RGCN_bases': lambda: RGCNConv(d, d, 967, num_bases=30).to(dev)(x, ei, et),
}
for name, fn in tests.items():
    try:
        out = fn()
        if out.requires_grad:
            out.sum().backward()
        torch.mps.synchronize(); print('OK  ', name)
    except Exception as e:
        print('FAIL', name, type(e).__name__, str(e)[:160])
print('recommended_max_memory GB:', torch.mps.recommended_max_memory() / 2**30)
```

### 6.2 Dense RGCN rewrite vs `RGCNConv` (CPU is enough)

```python
# rgcn_dense_check.py
import torch
from torch_geometric.nn import RGCNConv
torch.manual_seed(0)
n, R, B, d, E = 30, 7, 3, 8, 400
x = torch.randn(n, d)
ei = torch.randint(0, n, (2, E)); et = torch.randint(0, R, (E,))
conv = RGCNConv(d, d, R, num_bases=B)          # aggr='mean', root_weight=True, bias=True
ref = conv(x, ei, et)

A = torch.zeros(R, n, n)                        # A[r, i, j]: row = target i, column = source j
A.index_put_((et, ei[1], ei[0]), torch.ones(E), accumulate=True)
A = A / A.sum(-1, keepdim=True).clamp(min=1)    # per-relation mean, zero rows stay zero
M = torch.einsum('rb,rij->bij', conv.comp, A)
out = torch.einsum('bij,jk,bkl->il', M, x, conv.weight) + x @ conv.root + conv.bias

g_ref = torch.autograd.grad(ref.square().sum(), conv.comp)[0]
g_out = torch.autograd.grad(out.square().sum(), conv.comp)[0]
print('forward max abs diff:', (out - ref).abs().max().item())
print('grad(comp) max abs diff:', (g_out - g_ref).abs().max().item())
```

Both differences should be at float32 rounding level (~1e-6). This checks the algebra on all
relations at once; the drug-only split (5.1, step 3) then only restricts A to the drug rows.

### 6.3 Synthetic-size timing (replaces the (est.) figures until the Decagon data is downloaded)

Random edges at the Decagon counts of 4.1 (drug–drug edges among the first 645 nodes). Structure
differs from Decagon, so treat the result as an order-of-magnitude check, not a benchmark.

```python
# time_mps.py
import time, torch
from torch_geometric.nn import SAGEConv, GATConv, RGCNConv
dev = torch.device('mps'); N, D, d = 19_726, 645, 64

def graph(n_dd):
    dd = torch.randint(0, D, (2, n_dd))
    ppi = torch.randint(D, N, (2, 1_431_224))
    tg = torch.stack([torch.randint(0, D, (18_690,)), torch.randint(D, N, (18_690,))])
    ei = torch.cat([dd, ppi, tg, tg.flip(0)], 1)
    et = torch.cat([torch.randint(0, 964, (n_dd,)), torch.full((1_431_224,), 964),
                    torch.full((18_690,), 965), torch.full((18_690,), 966)])
    return ei.to(dev), et.to(dev)

def run(name, model, *args):
    x = torch.randn(N, d, device=dev)
    torch.mps.synchronize(); t = time.time()
    out = model(x, *args); out.square().mean().backward()
    torch.mps.synchronize()
    print(f'{name}: {time.time() - t:.2f} s fwd+bwd (1 layer), '
          f'driver memory {torch.mps.driver_allocated_memory() / 2**30:.1f} GB')
    torch.mps.empty_cache()

ei_c, _ = graph(71_000)
run('SAGEConv collapsed', SAGEConv(d, d).to(dev), ei_c)
run('GATConv collapsed', GATConv(d, d // 4, heads=4).to(dev), ei_c)
ei, et = graph(5_150_000)
run('RGCNConv 967 rel, bases=30', RGCNConv(d, d, 967, num_bases=30).to(dev), ei, et)
```

Run each encoder a second time and report the second timing (the first includes kernel
compilation). Expected [derived]: RGCN ≈ 5 GB for one layer; SAGE and GAT well under 1 GB.

### 6.4 Results on this Mac (run 2026-10-28)

Apple M4, 24 GB RAM (`torch.mps.recommended_max_memory()` = 16.0 GB), in the repo's `.venv`:
**torch 2.15.0** (not the pinned 2.14.1), torch_geometric 2.8.1. The scripts are the ones above,
unchanged except that 6.3 runs each encoder twice.

- **6.1 smoke test** [verified: run]: `scatter_reduce_amax`, `scatter_sum`, `scatter_max`,
  `softmax_index`, `SAGEConv`, `GATConv`, `HGTConv` and `RGCN_bases` all **OK** on MPS with no CPU
  fallback. `softmax_ptr` **FAIL**s with `NotImplementedError: The operator 'aten::segment_reduce'
  is not currently implemented for the MPS device`, exactly as predicted from the source.
- **6.2 dense rewrite** [verified: run]: forward max abs diff 9.5e-07; grad(comp) max abs diff
  4.6e-05. The rewrite matches `RGCNConv` at float32 rounding level.
- **6.3 synthetic-size timing** [measured, random edges at Decagon sizes, d = 64, one layer
  fwd + bwd, second run]:

  | Encoder | Time | MPS driver memory |
  |---|---|---|
  | SAGEConv, collapsed | 0.13 s | 1.0 GB |
  | GATConv, collapsed | 0.43 s | 3.0 GB |
  | RGCNConv, 967 relations, 30 bases | 3.28 s | 6.2 GB |
  | HGTConv, collapsed (ppi, targets, rev_targets, combo), heads = 4 | 0.53 s | 3.1 GB |

  HGT was timed separately (same conventions: random edges at Decagon sizes, 71,000 drug–drug
  edges, d = 64, second run; process peak 3.7 GB). It matches the derived ≈ 3–3.5 GB for the
  collapsed graph (4.5). HGT with all 967 edge types (derived ≈ 12–15 GB, 3,868-iteration Python
  loop without `pyg-lib`) was not run; 4.5 already rules it out.

  Driver memory includes the allocator's cache, so it is an upper bound on what the layer needs.
  It corrects two derived figures in section 4: **GAT measured 3.0 GB, not the ≈ 1–1.5 GB derived**
  (the first-pass 3–4 GB was closer), and **RGCN measured 6.2 GB vs. ≈ 5 GB derived**. Process peak
  memory footprint over the whole run: 7.3 GB. All four encoders fit comfortably in 16 GB.
- Still not measured: real epochs on the Decagon graph (needs the ppi / targets / mono files and a
  training script). Moved to the encoder comparison ticket #8 when #3 was closed.

## 7. What remains unverified

- **Runtime behaviour on this M4.** The `scatter_reduce_(amax)` MPS kernel is confirmed in the
  PyTorch source (2.6.0+, Metal from 2.13.0), but no op has been run here. Section 6.1 confirms it.
- **The dense rewrite.** The algebra matches the `RGCNConv` source (5.1), but it has not been run
  against `RGCNConv`. Section 6.2 does that.
- **All memory and time figures** in sections 1, 4 and 5 are derived from the source, not measured.
  Section 6.3 (synthetic sizes) and one timed epoch per encoder on the real data replace them.
- **Hardware figures:** M4 bandwidth (120 GB/s) and effective FP32 throughput (1 TFLOPS) are
  assumptions.
- Whether the latest `torch-scatter` / `torch-sparse` / `pyg-lib` releases accept MPS tensors. The
  evidence above (CUDA-only branches in PyG [S3], the maintainer statement [S9], wheels being CPU
  builds [S13]) says no.
- PyTorch 2.15.0 (on PyPI since 2026-10-28) was not checked against its source (see Blockers).
- The MPS backend's float64 support (5.5).

## 8. Blockers

- **Resolved 2026-10-28 (see 6.4):** the user set up `.venv` in the repo; 6.1–6.3 were run.
  Original blocker: **No PyTorch / PyG on this Mac.** `/Library/Frameworks/Python.framework/Versions/3.13/bin/python3`
  has neither `torch` nor `torch_geometric`, and nothing was installed in this research. The smoke
  test (6.1), the rewrite check (6.2) and the timing (6.3) were therefore not run. **User action:**
  run the install commands in section 6, then `python smoke_mps.py`, `python rgcn_dense_check.py`
  and `python time_mps.py`, and paste the output into issue #3.
- **Decagon data not downloaded.** Real timed epochs per encoder are not possible yet. **User
  action:** after the data download (`3.1`), time one epoch per encoder with the project's training
  script.
- **PyTorch 2.15.0 source not readable.** PyPI lists `torch` 2.15.0 (uploaded
  2026-10-28T15:16:30), but
  `curl https://raw.githubusercontent.com/pytorch/pytorch/v2.15.0/aten/src/ATen/native/mps/operations/ScatterGather.mm`
  returned **HTTP 404** (no `v2.15.0` tag yet), and the GitHub releases list ends at v2.14.1. So
  the install above pins `torch==2.14.1`.
- **M4 memory bandwidth not re-confirmed.** `https://www.apple.com/macbook-air/specs/` and
  `https://www.apple.com/macbook-pro/specs/` (HTTP 200) now list M5 models (153 GB/s and up);
  `https://www.apple.com/newsroom/2024/05/apple-introduces-m4-chip/` (HTTP 200) gives no GB/s
  figure. The 120 GB/s used in the time estimates stays unverified.

## Sources

- **[S1]** PyTorch 2.14 docs, "MPS backend" (macOS 14.0+ requirement). https://docs.pytorch.org/docs/2.14/notes/mps.html
- **[S2]** PyG GitHub releases (2.8.1 latest; 2.8.0: torch-cluster/spline-conv replaced by `pyg-lib==0.7.0`; 2.6.0: EdgeIndex to deprecate `SparseTensor`). https://github.com/pyg-team/pytorch_geometric/releases
- **[S3]** PyG source, `torch_geometric/utils/_scatter.py` @ `ea2c5f6`. https://github.com/pyg-team/pytorch_geometric/blob/ea2c5f62d54c8b414dd3ea7139814a7c1ef83182/torch_geometric/utils/_scatter.py#L62-L138
- **[S4]** PyG source, `torch_geometric/utils/_softmax.py` @ `ea2c5f6` (ptr path L55–81, index path with max from `src.detach()` L82–88). https://github.com/pyg-team/pytorch_geometric/blob/ea2c5f62d54c8b414dd3ea7139814a7c1ef83182/torch_geometric/utils/_softmax.py#L55-L92
- **[S5]** PyTorch 2.14 docs, "MPS Environment Variables". https://docs.pytorch.org/docs/2.14/mps_environment_variables.html
- **[S6]** PyG source, `torch_geometric/nn/dense/linear.py` @ `ea2c5f6` (`HeteroLinear.forward_naive` L248–252, `forward` L287). https://github.com/pyg-team/pytorch_geometric/blob/ea2c5f62d54c8b414dd3ea7139814a7c1ef83182/torch_geometric/nn/dense/linear.py#L248-L252
- **[S7]** PyG source and docstring, `torch_geometric/nn/conv/rgcn_conv.py` @ `ea2c5f6` (`masked_edge_index` L23–26, defaults L99–100, bases L205, `segment_matmul` condition L242–246, relation loop L257–270, root/bias L275–280, `FastRGCNConv` L314–374). https://github.com/pyg-team/pytorch_geometric/blob/ea2c5f62d54c8b414dd3ea7139814a7c1ef83182/torch_geometric/nn/conv/rgcn_conv.py
- **[S8]** PyG installation docs (optional extensions, wheel matrix, no conda after PyTorch 2.5; `pip install torch_geometric`; Python 3.10–3.14). https://pytorch-geometric.readthedocs.io/en/latest/install/installation.html and https://github.com/pyg-team/pytorch_geometric/blob/ea2c5f62d54c8b414dd3ea7139814a7c1ef83182/docs/source/install/installation.rst
- **[S9]** PyG issue 4776, "Is it possible to get GPU acceleration support on the M1 chip?" (rusty1s, 2022-06-09). https://github.com/pyg-team/pytorch_geometric/issues/4776
- **[S10]** PyTorch 2.13 release blog, "Large MPS Op Migration to Native Metal". https://pytorch.org/blog/pytorch-2-13-release-blog/
- **[S11]** PyG issue 6837, "Training TGN with MPS not supported on Mac M1 max" (`aten::_unique2`). https://github.com/pyg-team/pytorch_geometric/issues/6837
- **[S12]** PyG issue 9754, "PYG compatible with Metal pytorch 2.5 version". https://github.com/pyg-team/pytorch_geometric/issues/9754
- **[S13]** PyG wheel index for PyTorch 2.8 (macOS wheel filenames). https://data.pyg.org/whl/torch-2.8.0+cpu.html
- **[S14]** PyTorch GitHub releases (v2.14.1 2026-09-30, v2.14.0 2026-09-02, v2.13.0 2026-07-08). https://github.com/pytorch/pytorch/releases
- **[S15]** PyTorch `aten/src/ATen/native/native_functions.yaml`, entry `scatter_reduce.two_out`:
  - [S15a] v2.5.0 L8250 `CPU, CUDA: scatter_reduce_two`. https://github.com/pytorch/pytorch/blob/v2.5.0/aten/src/ATen/native/native_functions.yaml#L8246-L8250
  - [S15c] v2.6.0 L8291 `CPU, CUDA, MPS: scatter_reduce_two`. https://github.com/pytorch/pytorch/blob/v2.6.0/aten/src/ATen/native/native_functions.yaml#L8287-L8291
  - v2.14.0 L8368 `CPU, CUDA, MPS, XPU: scatter_reduce_two`. https://github.com/pytorch/pytorch/blob/v2.14.0/aten/src/ATen/native/native_functions.yaml#L8364-L8368
  - [S15b] v2.14.0 L14868–14872, `segment_reduce` dispatch `CPU, CUDA, XPU` (no MPS). https://github.com/pytorch/pytorch/blob/v2.14.0/aten/src/ATen/native/native_functions.yaml#L14868-L14872
  - [S15d] v2.14.0 `aten/src/ATen/native/TensorAdvancedIndexing.cpp` (CPU-only expanded-index path L2037, `scatter_reduce_exclude_self_helper` L2133–2160, `scatter_impl` L2232–2263, `scatter_reduce_two` L2381–2425). https://github.com/pytorch/pytorch/blob/v2.14.0/aten/src/ATen/native/TensorAdvancedIndexing.cpp#L2232-L2263
- **[S16]** PyTorch PR #141948 "[MPS] Add scatter_reduce.two" (malfet; commit `d6481333ad33`, 2024-12-04; in v2.6.0, not in v2.5.0). https://github.com/pytorch/pytorch/pull/141948
- **[S17]** PyTorch PR #184028 "[MPS] Migrate scatter/gather ops to Metal with async bounds checking" (Nikita Shulga; commit `9b7a45e57afd`, 2026-05-18; in v2.13.0, not in v2.12.0). https://github.com/pytorch/pytorch/pull/184028
  - [S17a] v2.6.0 `aten/src/ATen/native/mps/operations/ScatterGather.mm` L343–367 (MPSGraph-based `scatter_reduce_two_mps_kernel`, `amax` case). https://github.com/pytorch/pytorch/blob/v2.6.0/aten/src/ATen/native/mps/operations/ScatterGather.mm#L343-L367
  - [S17b] v2.14.0 `ScatterGather.mm` (`scatter_reduce_metal` L189, `scatter_reduce_dispatch` L439–473, `scatter_reduce_two_mps_kernel` + `REGISTER_DISPATCH` L488–501). https://github.com/pytorch/pytorch/blob/v2.14.0/aten/src/ATen/native/mps/operations/ScatterGather.mm#L439-L501
  - [S17c] v2.14.0 `aten/src/ATen/native/mps/kernels/ScatterGather.metal` (amax via `atomic_binary_op` L486–495; `REGISTER_SCATTER_OPS_FULL(float)` L519). https://github.com/pytorch/pytorch/blob/v2.14.0/aten/src/ATen/native/mps/kernels/ScatterGather.metal#L481-L530
- **[S18]** PyTorch v2.14.0 `torch/testing/_internal/common_mps.py`, `UNIMPLEMENTED_XFAILLIST` (L368) with `_segment_reduce*` (L586–587); only `scatter_reduceprod` int64 listed for scatter_reduce (L585). https://github.com/pytorch/pytorch/blob/v2.14.0/torch/testing/_internal/common_mps.py#L580-L590
- **[S19]** PyG source, `torch_geometric/nn/conv/hgt_conv.py` @ `ea2c5f6` (`num_types` L79–83, per-type source features L119–157, bipartite propagate L195–199, softmax L230). https://github.com/pyg-team/pytorch_geometric/blob/ea2c5f62d54c8b414dd3ea7139814a7c1ef83182/torch_geometric/nn/conv/hgt_conv.py
- **[S20]** PyG source, `torch_geometric/nn/conv/gat_conv.py` @ `ea2c5f6` (softmax L404, message L408–409). https://github.com/pyg-team/pytorch_geometric/blob/ea2c5f62d54c8b414dd3ea7139814a7c1ef83182/torch_geometric/nn/conv/gat_conv.py#L387-L409
- **[S21]** PyG source, `torch_geometric/nn/conv/message_passing.py` @ `ea2c5f6` (`ptr` handling L379–400). https://github.com/pyg-team/pytorch_geometric/blob/ea2c5f62d54c8b414dd3ea7139814a7c1ef83182/torch_geometric/nn/conv/message_passing.py#L379-L400
- **[S22]** PyG source, `torch_geometric/utils/_segment.py` @ `ea2c5f6` (`segment` L27–34, `_torch_segment` → `torch._segment_reduce` L37–50). https://github.com/pyg-team/pytorch_geometric/blob/ea2c5f62d54c8b414dd3ea7139814a7c1ef83182/torch_geometric/utils/_segment.py#L27-L50
- **[S23]** PyTorch v2.14.0 `torch/mps/__init__.py` (`synchronize` L30, `empty_cache` L78, `current_allocated_memory` L110, `driver_allocated_memory` L120, `recommended_max_memory` L130). https://github.com/pytorch/pytorch/blob/v2.14.0/torch/mps/__init__.py
- **[S24]** PyG source, `torch_geometric/nn/aggr/base.py` (`Aggregation.reduce` L172–185) and `torch_geometric/nn/aggr/basic.py` (`MeanAggregation.forward` → `reduce='mean'`) @ `ea2c5f6`. https://github.com/pyg-team/pytorch_geometric/blob/ea2c5f62d54c8b414dd3ea7139814a7c1ef83182/torch_geometric/nn/aggr/base.py#L172-L185
