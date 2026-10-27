# PyG on Apple Silicon (MPS): feasibility of the Decagon encoders

Research for issue #3 (map #1). Researched 2026-10-27.
Point IDs (`3.2`, `5.3`, ...) refer to `PROJECT_CONTEXT.md`.

**Machine checked:** Apple M4, 24 GB unified memory, macOS 15.5 (`sysctl`, `sw_vers`).
**Software state at time of writing:** PyTorch 2.14 docs are current [S1]; PyG 2.8.1 is the latest release [S2].

Legend: **[Sn]** = cited primary source (list at the end). **(est.)** = my own arithmetic from the
Decagon counts in `3.1`/`3.2`, not a measurement. **(unverified)** = could not be confirmed from a
primary source during this research; settle it with the smoke test in section 6.

---

## 1. Short answer

| Encoder (`5.3`) | Runs on MPS? | Memory, 2 layers, d = 64 (est.) | Locally feasible on the M4 / 24 GB? |
|---|---|---|---|
| Tensor factorisation (no encoder) | Yes. Only `nn.Embedding`, elementwise ops and matmul. | < 0.5 GB | **Yes**, easily. Seconds per epoch. |
| `to_hetero` GraphSAGE, collapsed drug–drug edges | Yes in principle: sum/mean aggregation uses plain `scatter_add_` / `index_select` [S3]. | ~1–3 GB | **Yes.** Seconds to tens of seconds per epoch. |
| `to_hetero` GAT, collapsed drug–drug edges | Needs a scatter **max** (softmax) → `torch.scatter_reduce_(reduce='amax')` on MPS [S3, S4]. MPS kernel status (unverified); CPU fallback exists [S5]. | ~3–4 GB | **Probably yes**; confirm with the smoke test. |
| HGT, collapsed drug–drug edges (4 edge types) | Same softmax dependency as GAT; type-specific linears fall back to a plain Python loop without `pyg-lib` [S6]. | ~3–5 GB | **Probably yes**, if collapsed. With 964 edge types: ~13 GB and a 964-step loop, not practical. |
| `RGCNConv`, 967 relations, basis decomposition | Runs with core ops, but loops over all 967 relations per layer; `segment_matmul` path is disabled when `num_bases` is set [S7]. | **~10 GB** (saved per-relation outputs) | **Marginal.** Fits only just; slow (several seconds per encoder pass, est.). Better with the rewrite in 5.1. |
| `FastRGCNConv`, same setup | Same ops, but materialises one weight matrix per edge [S7]. | **~108 GB** per layer | **No.** Infeasible on any 24 GB machine. |

Bottom line: tensor factorisation, GraphSAGE and (likely) GAT/HGT on the collapsed graph run locally
on MPS. Stock `RGCNConv` with 967 relations is borderline; `FastRGCNConv` is not feasible. A
mathematically equivalent rewrite of the RGCN drug–drug block (5.1) makes the reference encoder cheap
and fully dense, which is the best-supported kind of work on MPS.

---

## 2. How PyG executes on MPS

- **Core PyG needs no extension libraries.** "From PyG 2.3 onwards, you can install and use PyG
  without any external library required except for PyTorch" [S8]. In that mode, message passing is
  built from core PyTorch ops, so it runs on whatever device the tensors are on, including MPS.
- **PyG's maintainer said this early on.** Core ops should get GPU acceleration on M1 for layers such
  as `GCNConv` / `SAGEConv` that do not rely on `torch-scatter` / `torch-sparse`. MPS support for
  PyG's own external operators is "of lower priority" [S9, rusty1s, 2022-06-09].
- **Which core op each aggregation uses** (PyG `utils/_scatter.py`, current master) [S3]:
  - `sum` / `add` / `mean`: `scatter_add_` (and `scatter_` for `any`). There are no device checks.
  - `min` / `max`: `torch.scatter_reduce_(reduce='amax'|'amin', include_self=False)` unless the
    input is **CUDA** and requires grad, and `torch-scatter` is installed. On MPS this means the
    core PyTorch kernel is always used.
  - `mul`: `scatter_reduce_(reduce='prod')`.
  - The only device branch in the file is `src.is_cuda`; there is **no MPS-specific code path** in
    PyG's scatter.
- **Attention layers need a scatter max.** PyG's `softmax` subtracts a per-group max taken from
  `src.detach()` before `exp()` [S4]. GAT and HGT therefore depend on MPS support for
  `scatter_reduce_(amax)`. Because the input is detached, this path is used even when
  `torch-scatter` is installed.
- **PyTorch MPS status.** The MPS backend needs macOS 14.0+ [S1]. PyTorch 2.13 moved "scatter/gather
  with bounds checking", sort, reductions and embedding backward from MPSGraph to hand-written Metal
  kernels [S10]. I could not find a primary source that confirms whether `scatter_reduce` with
  `amax` is natively implemented on MPS (unverified).
- **CPU fallback for missing ops.** `PYTORCH_ENABLE_MPS_FALLBACK=1` makes PyTorch "fall back to the
  CPU for operations that MPS doesn't support" [S5]. Without it, a missing op raises
  `NotImplementedError ... is not currently implemented for the MPS device`. An example is
  `aten::_unique2` in PyG's TGN example [S11].
- **Memory limits on MPS** [S5]:
  - `PYTORCH_MPS_HIGH_WATERMARK_RATIO` (default 1.7) is a hard cap, as a multiple of the device's
    recommended maximum working set.
  - `PYTORCH_MPS_LOW_WATERMARK_RATIO` (default 1.4 on unified memory) is a soft limit.
  - On a 24 GB machine the working set is shared with the CPU side; query it with
    `torch.mps.recommended_max_memory()`.
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
- Supervision per epoch: ≈ 1.10M positives + 1.10M degree-matched negatives = **2.2M scored
  triples**.
- Tensor factorisation has no message graph, so all training triples are supervision:
  ≈ 7.4M per epoch.

**Model:**
- d = 64, float32, 2 encoder layers.
- Decagon decoder costs ≈ 2d² + 4d ≈ 8.4 kFLOP per triple forward, ≈ 25 kFLOP forward + backward.

**Hardware:**
- M4 memory bandwidth is 120 GB/s (Apple spec).
- FP32 GPU throughput: Apple does not publish it; ~4 TFLOPS peak is a commonly quoted figure.
- I assume 1 TFLOPS effective in the time estimates.

All numbers below are (est.).

### 4.2 Tensor factorisation (no encoder)

- **Ops:** `Embedding`, `index_select`, elementwise, `matmul`, `BCEWithLogitsLoss`. All are standard
  MPS ops.
- **Parameters:** 645·64 + 964·64 + 64² ≈ 0.11M.
- **Time:** ≈ 7.4M × 25 kFLOP ≈ 0.19 TFLOP per epoch, so well under a second of compute. In practice
  per-batch overhead dominates: with ~64k-triple batches that is about 115 steps, so a few seconds
  per epoch.
- **Memory:** ≈ 1.5 KB of activations per triple, so ≈ 100 MB per 64k batch.
- **Verdict:** Feasible on MPS; the cheapest model.

### 4.3 GraphSAGE via `to_hetero`, drug–drug edges collapsed to one type

- **Ops:** gather (`index_select`), `scatter_add_` / mean, `Linear`. All core ops [S3].
- **Edges per layer:** 1.43M PPI + 37k target + 71k drug–drug ≈ 1.54M. PPI dominates; the
  drug–drug part is tiny once collapsed.
- **Memory:** about 400 MB of transient messages per layer, so ~1–3 GB peak for 2 layers
  forward + backward.
- **Time:** a single encoder pass moves ~1–2 GB of memory, ~0.1–0.5 s with overhead.
- **Per epoch, two ways to train:**
  - Full-batch: one encoder pass plus the decoder over 2.2M triples in chunks, then one backward.
    About 1–3 s, with ~3 GB of decoder activations.
  - Minibatch: re-run the encoder per 64k batch, about 34 passes, so ~5–20 s.
- **Verdict:** Feasible.

### 4.4 GAT via `to_hetero`, collapsed

- **Ops:** the same as GraphSAGE, plus `softmax`, which needs a scatter max
  (`scatter_reduce_(amax)`) [S3, S4].
- **MPS status:** if that kernel is missing on MPS it raises `NotImplementedError`, or with
  `PYTORCH_ENABLE_MPS_FALLBACK=1` it runs on the CPU per call (unverified).
- **Memory:** roughly 1 KB saved per edge per layer (messages, attention coefficients, dropout
  mask), so ~1.5 GB per layer and ~3–4 GB in total.
- **Time:** about 2–3× GraphSAGE.
- **Practical caveat:** `GATConv` on bipartite edge types (drug→protein) needs
  `add_self_loops=False` under `to_hetero`.
- **Verdict:** Likely feasible; one smoke test settles it.

### 4.5 HGT, collapsed (edge types: ppi, targets, rev_targets, combo)

- **Ops:**
  - Q/K/V projections per node type and relation transforms per edge type.
  - Attention softmax (scatter max) as in GAT.
- **Without `pyg-lib`:**
  - `HeteroLinear` falls back to `forward_naive`, a Python loop over types doing
    `x[start:end] @ weight[i]` [S6].
  - `HeteroDictLinear` loops over node types with ordinary `Linear` modules [S6].
  - Both are fine on MPS for 2 node types and 4 edge types.
- **Memory:** ~1 KB saved per edge per layer, so ~3–5 GB.
- **With all 964 side-effect edge types kept:**
  - E ≈ 6.6M edges, about 13 GB for 2 layers.
  - A ~967-iteration type loop with a data-dependent slice in each iteration.
  - Not practical locally.
- **Verdict:** Feasible only with collapsed edge types; depends on the same softmax check as GAT.

### 4.6 `RGCNConv`, 967 relations, `num_bases = 30`, on `data.to_homogeneous()`

How the PyG layer runs in this configuration [S7]:

- **Bases.** The basis weights are expanded once per forward: `comp @ basis` gives
  [967, 64, 64] ≈ 16 MB.
- **`segment_matmul` is off.** It is used only if `pyg-lib` provides it, the layer is not being
  compiled, **`num_bases` is unset**, and the input is floating point. With bases it is always off.
  On MPS it would be unavailable anyway.
- **Loop over relations.** The layer runs a **Python loop over all 967 relations**. In each
  iteration it:
  - builds a boolean mask (`edge_type == i`) over all edges and applies it;
  - runs `propagate` on that relation's edges;
  - computes `out = out + h @ weight[i]`.
- **The docs warn about this pattern.** `RGCNConv` "can use dynamic shapes" and "may result in low
  GPU utilization in case the graph has a large number of relations".

Cost (est.), with E ≈ 6.62M edges in the homogeneous graph:

- **Masking.**
  - Each iteration reads the full `edge_index` (≈ 106 MB), so about 100 GB of memory traffic per
    layer forward. That is ≈ 0.9 s per layer at 120 GB/s.
  - Each iteration also has a data-dependent output size: 967 device synchronisations per layer.
- **Saved activations.**
  - For the backward of `h @ weight[i]`, every relation saves `h` of shape [N, d] with
    N = 19,726, even though 964 of the relations only touch the 645 drugs.
  - 967 × 19,726 × 64 × 4 B ≈ **4.9 GB per layer, ≈ 10 GB for 2 layers** (d = 32: ≈ 5 GB).
- **Time.** About 2–10 s per encoder forward + backward. Per epoch that is a few seconds to tens of
  seconds full-batch, or minutes with per-batch encoder passes.
- **Verdict:** Runs on MPS with core ops, but uses about 40–60 % of the realistic MPS budget on
  24 GB and is slow. Marginal.

### 4.7 `FastRGCNConv`, same setup

- **How it works:** one `propagate` over all edges with no loop over relations. Its `message` does a
  per-edge batched matmul with `weight[edge_type]` [S7]. The docs warn it "may consume a large
  amount of memory to compensate" [S7].
- **Memory (est.):**
  - The gathered weight tensor is [E, 64, 64]: 6.62M × 64 × 64 × 4 B ≈ **108 GB per layer**.
  - Drug–drug edges alone: ≈ 84 GB.
  - Even at d = 16: ≈ 7 GB per layer.
- **Verdict:** Infeasible locally, and on most single GPUs.

## 5. Workarounds

### 5.1 Rewrite RGCN's drug–drug block as dense drug-only algebra (derivation, not from a source)

The 964 side-effect relations connect only the 645 drugs. With basis decomposition
W_r = Σ_b a_rb V_b, the relational sum for drug i is a reordering of the same sum:

Σ_r Σ_{j ∈ N_r(i)} (1/c_ir) W_r x_j = Σ_b [ (Σ_r a_rb Â_r) X_drug ]_i V_b

Here Â_r is the row-normalised 645 × 645 adjacency of side effect r, so the per-relation mean
normalisation is folded in.

**Implementation:**
1. Precompute `A = stack(Â_r)` as a constant dense tensor [964, 645, 645]. This is ≈ 1.6 GB in
   float32, or 0.8 GB in float16.
2. Each layer:
   - `M = einsum('rb,rij->bij', comp, A)` gives [30, 645, 645]. That is ≈ 24 GFLOP forward and
     about as much again for the gradient with respect to `comp`.
   - `out_drug = Σ_b M_b @ X_drug @ V_b`.
3. Handle the 3 protein relations (PPI, targets, reverse targets) with ordinary PyG layers, or
   `RGCNConv(num_relations=3)`, and add the results.

**Properties:**
- Mathematically identical to `RGCNConv` with bases on these relations.
- Uses only dense `einsum`/`matmul`, the best-supported work on MPS.
- ≈ 2 GB of memory and tens of milliseconds per layer (est.).
- Removes the 967-iteration loop and the ~10 GB of saved activations.
- **Cost:** it is a custom layer, so it needs its own test against `RGCNConv` on a small graph.

### 5.2 Message passing on the CPU, decoder (and TF) on MPS

- **How:** autograd tracks `.to()` across devices. `z = encoder_cpu(...)` then `z.to('mps')`, with
  the loss on MPS and `loss.backward()`, works.
- **Transfer cost:** only z_drug [645, 64] ≈ 165 KB per step.
- **Why it can help:**
  - PyG's own comment says `torch.scatter_reduce` is faster on CPU than `torch-scatter` [S3].
  - On the CPU, `pyg-lib` CPU kernels (e.g. `segment_matmul` for `HeteroLinear`) become usable [S6].
- **Costs:**
  - The CPU and GPU share the same 120 GB/s unified memory bandwidth.
  - Bandwidth-bound scatter work may therefore not speed up much in either place.

### 5.3 `PYTORCH_ENABLE_MPS_FALLBACK=1`

- **Effect:** any op MPS lacks (e.g. `aten::_unique2` [S11]) runs on the CPU instead of raising an
  error [S5]. Slower, but correct.
- **How to set it:** in the environment before starting Python.
- **Debugging tip:** leave it **off** while testing, so missing ops show up as errors.

### 5.4 Subsample side effects

- **Effect on RGCN:** loop length, saved activations and masked-edge traffic all scale linearly with
  the number of relations. Top 100 side effects: ≈ 0.5 GB per layer instead of 4.9 GB (est.).
- **Cost:** results are then no longer comparable with Decagon's 964-type numbers (`6.2`). Use it
  for development and debugging, not for the final comparison.

### 5.5 Smaller d, full-batch training, mixed precision

- **Smaller d:** d = 32 halves all activation figures above.
- **Full-batch training:** one encoder pass per epoch, with the decoder loss accumulated over chunks
  of supervision triples. This removes the dominant "encoder re-run per batch" factor.
- **Mixed precision:** float16 for the constant `A` in 5.1 halves its memory. float64 should be
  avoided entirely. A float64 multi-hot drug feature matrix from pandas must be cast to float32
  before `.to('mps')`. The MPS backend does not support float64 (unverified: not checked against a
  primary source here).

## 6. Smoke test to settle the unverified points

Run once on the M4 with **no** `PYTORCH_ENABLE_MPS_FALLBACK` set. Any `NotImplementedError` names the
missing op.

```python
import torch
from torch_geometric.nn import SAGEConv, GATConv, HGTConv, RGCNConv
from torch_geometric.utils import scatter, softmax

dev = torch.device('mps')
N, E, d = 2000, 20000, 64
x = torch.randn(N, d, device=dev)
ei = torch.randint(0, N, (2, E), device=dev)
et = torch.randint(0, 967, (E,), device=dev)
idx = ei[1]

for name, fn in {
    'scatter_sum': lambda: scatter(x[ei[0]], idx, 0, dim_size=N, reduce='sum'),
    'scatter_max': lambda: scatter(x[ei[0]], idx, 0, dim_size=N, reduce='max'),
    'softmax':     lambda: softmax(torch.randn(E, device=dev), idx, num_nodes=N),
    'SAGEConv':    lambda: SAGEConv(d, d).to(dev)(x, ei),
    'GATConv':     lambda: GATConv(d, d // 4, heads=4).to(dev)(x, ei),
    'RGCN_bases':  lambda: RGCNConv(d, d, 967, num_bases=30).to(dev)(x, ei, et),
}.items():
    try:
        out = fn(); out.sum().backward() if out.requires_grad else None
        torch.mps.synchronize(); print('OK  ', name)
    except Exception as e:
        print('FAIL', name, type(e).__name__, str(e)[:120])
print('recommended_max_memory GB:', torch.mps.recommended_max_memory() / 2**30)
```

Add an `HGTConv` check on a two-type `HeteroData` if HGT stays in the comparison.

## 7. What remains unverified

- Whether `scatter_reduce_(amax)` (GAT/HGT softmax) has a native MPS kernel in PyTorch 2.13/2.14,
  as opposed to needing the CPU fallback. The 2.13 notes confirm native Metal scatter/gather [S10],
  but not the reduce variants.
- Whether the latest `torch-scatter` / `torch-sparse` / `pyg-lib` releases accept MPS tensors. The
  evidence above (CUDA-only branches in PyG [S3], the maintainer statement [S9], wheels being CPU
  builds [S13]) says no.
- Real wall-clock times. All times above are order-of-magnitude estimates; the smoke test plus one
  timed epoch replaces them.

## Sources

- **[S1]** PyTorch 2.14 docs, "MPS backend" (macOS 14.0+ requirement). https://docs.pytorch.org/docs/2.14/notes/mps.html
- **[S2]** PyG GitHub releases (2.8.1 latest; 2.8.0: torch-cluster/spline-conv replaced by `pyg-lib==0.7.0`; 2.6.0: EdgeIndex to deprecate `SparseTensor`). https://github.com/pyg-team/pytorch_geometric/releases
- **[S3]** PyG source, `torch_geometric/utils/_scatter.py` (master). https://github.com/pyg-team/pytorch_geometric/blob/master/torch_geometric/utils/_scatter.py
- **[S4]** PyG source, `torch_geometric/utils/softmax` (max from `src.detach()`). https://pytorch-geometric.readthedocs.io/en/2.3.0/_modules/torch_geometric/utils/softmax.html
- **[S5]** PyTorch 2.14 docs, "MPS Environment Variables". https://docs.pytorch.org/docs/2.14/mps_environment_variables.html
- **[S6]** PyG source, `torch_geometric/nn/dense/linear.py` (`HeteroLinear.forward_naive`, `HeteroDictLinear`). https://github.com/pyg-team/pytorch_geometric/blob/master/torch_geometric/nn/dense/linear.py
- **[S7]** PyG source and docstring, `torch_geometric/nn/conv/rgcn_conv.py` (`RGCNConv`, `FastRGCNConv`). https://github.com/pyg-team/pytorch_geometric/blob/master/torch_geometric/nn/conv/rgcn_conv.py
- **[S8]** PyG installation docs (optional extensions, wheel matrix, no conda after PyTorch 2.5). https://pytorch-geometric.readthedocs.io/en/latest/install/installation.html
- **[S9]** PyG issue 4776, "Is it possible to get GPU acceleration support on the M1 chip?" (rusty1s, 2022-06-09). https://github.com/pyg-team/pytorch_geometric/issues/4776
- **[S10]** PyTorch 2.13 release blog, "Large MPS Op Migration to Native Metal". https://pytorch.org/blog/pytorch-2-13-release-blog/
- **[S11]** PyG issue 6837, "Training TGN with MPS not supported on Mac M1 max" (`aten::_unique2`). https://github.com/pyg-team/pytorch_geometric/issues/6837
- **[S12]** PyG issue 9754, "PYG compatible with Metal pytorch 2.5 version". https://github.com/pyg-team/pytorch_geometric/issues/9754
- **[S13]** PyG wheel index for PyTorch 2.8 (macOS wheel filenames). https://data.pyg.org/whl/torch-2.8.0+cpu.html
