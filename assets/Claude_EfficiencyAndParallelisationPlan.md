# Step-matrix efficiency and parallel row computation

## Context
Building the step matrix S dominates `PathIntegration` run time. For example, d=2 Duffing with QuinticAxis, N=81, RK4 takes 1.06 s serially.

A profile shows where the time goes:
- ~55% in `temp .*= w; res .+= temp` over the full N^d array for every quadrature node (discreteintegrator.jl:83-86);
- ~23% in Newton back-tracing;
- ~5% in sparse `setindex!`, which shifts `colptr` on every insertion.

Sparse interpolations therefore cost O(N^d) per row and O(N^{2d}) overall. That is an implementation artefact: the method only needs O(n_q·(o+1)^d) per row.

Each row `Φ_(i)` (paper Eq. 21–22) is independent of the others, so rows can also be computed in parallel. That needs per-task caches and a race-free sparse assembly.

This plan does the efficiency work first (Part A), then the parallelisation (Part B), which reuses A1–A3. Everything is in the serial path too.

---

## Part A — Efficiency improvements (ranked)

### A1. Sparse accumulator (SPA) row kernel for sparse interpolations
**Workspace.** The per-task value buffer is the existing `IK.temp.itpM` (N^d, all zero between rows). Add a `mark::Vector{Bool}` (N^d) and `touched::Vector{Int}` for the indices touched in the current row.

**Per row:**
1. Reset only the previous row's `touched` entries (`itpM[j] = 0`, `mark[j] = false`, then `empty!`).
2. If `Q_integrate[]`, loop over the nodes `(w, x)`:
   - `fx = kernel_value!(IK, x)` (TPDF·detJ, or 0 below 1e-8);
   - `basefun_vals_safe!(IK)`;
   - for each stencil entry: `v = (fx*prod(val))*w`, skip if `iszero(v)`, compute the linear index `j`, mark it if new, and do `itpM[j] += v`.
3. The row output does `sort!(touched)` and emits entries in ascending j.

**Cost:** O(n_q·(o+1)^d + nnz_row·log nnz_row) per row, instead of O(n_q·N^d).

**Memory:** `di.temp` and `di.res` (N^d each) are no longer needed and are created with size 0. Per task this goes from 24 B·N^d to 9 B·N^d.

**Exactness:** the arithmetic matches the old path: `(fx*prod(val))*w` is what `temp .*= w` computes, nodes are summed in the same order, and adding +0.0 is a no-op. So SPA rows are **bitwise equal** to the old generic kernel. A test checks this.

**Where it applies:** pdfs that `get_itp_type` calls sparse (all-sparse or mixed axes) together with a `DiscreteIntegrator{1}`. QuadGK and the vibro-impact `NonSmoothDiscreteIntegrator` keep the generic kernel.

**Code:**
- Split the kernel functor `(IK)(vals,x)` (integrationkernel.jl:17-32) into `kernel_value!(IK,x)` plus the existing basis/fill steps. The generic functor calls it, so its behaviour is unchanged.
- Add a `kernel` field to `IK_temp`, holding one of:
  - `GenericRowKernel()`
  - `SparseAccumulator(mark, touched)`
  - `DenseTensorKernel(…)` (A4)
- `get_IK_weights!` dispatches on the kernel type.
- The kernel is chosen in `PathIntegration` (pathintegration.jl:77-87), which must now build `di` before `ikt`.
- An internal keyword `generic_row_kernel = false` (undocumented) forces the old kernel, for tests and benchmarks.

### A2. Direct CSC assembly (plus optional `Int32` indices)
Rows of S are columns of the stored Sᵀ, so the matrix is assembled directly:
- For each row, `fill_to_stepMX!` appends `(rowval, nzval)` to a row buffer in ascending order and records `nnz_per_row[i]`. With SPA it iterates the sorted `touched`; otherwise it scans `itpM`.
- Then `colptr = 1 .+ cumsum`, followed by `resize!`/`copyto!` into the existing CSC. This unwraps `Transpose` → `parent` and `ThreadedSparseMatrixCSC` → `.A`, so `recompute_stepMX!` fills the wrapped matrix in place.

This is O(nnz): it removes the quadratic `setindex!` insertion, and it is the same code the chunked parallel version needs (Part B).

`SparseMX(; index_type = Int)`:
- Passing `Int32` stores 12 B instead of 16 B per nonzero, which means less memory traffic per time step.
- `initialize_stepMX` uses `spzeros(T, Ti, l, l)`.
- Assembly throws an informative `ArgumentError` if nnz would overflow `Ti`.
- The default stays `Int`, because d=4 can exceed 2³¹ nonzeros.

### A3. Stateless quadrature mapping
`DiscreteIntegrator` gains `x_ref::xT, w_ref::wT`. `rescale_xw!` maps `[x_ref[1], x_ref[end]] → [start, stop]` from the reference nodes every time, instead of from the previous row's nodes. This:
- removes the accumulated rescaling error;
- makes each row depend only on its grid point, which is needed for thread-count-independent, bitwise results;
- costs nothing.

The constructor also `collect`s `x`, so `Trapezoidal`/`NewtonCotesIntegrator` (`x::LinRange`) no longer fail under `smart_integration` (a latent bug). Files: `src/types.jl`, `src/integration/discreteintegrator.jl:13-20,138-166`.

### A4. Dense interpolations: vectorised row
**Workspace:** `DenseTensorKernel(Bs, c, B1c, KR)`, where each `Bⱼ` is Nⱼ×n_q.

**Per node q:**
- `c[q] = w_q·fx_q`, or 0 with a zeroed column if fx is negligible;
- `basefun_vals_safe!(view(Bⱼ,:,q), axisⱼ, x0[j])`. The Chebyshev and trigonometric basis functions already work on views.

**Then one BLAS product:**
- d=1: `mul!(vec(itpM), B₁, c)`;
- d=2: `B1c .= B₁ .* c'` then `mul!(itpM, B1c, transpose(B₂))` (Φ = B₁·diag(c)·B₂ᵀ, gemm);
- d≥3: build the Khatri–Rao matrix `KR[:,q] = kron(B_d[:,q],…,B₂[:,q])`, then `mul!(reshape(itpM,N₁,:), B1c, transpose(KR))`.

This applies when all axes are dense and the integrator is a `DiscreteIntegrator{1}`. It changes the summation order, so results match the generic kernel to about 1e-13 relative.

### A5. Direct steady-state solver
New exported `steady_state!(PI; method = :auto, tol = 1e-10, reset_t = true) -> (PI, info)`.

**Operator:** S, or for time-periodic systems the period map P = S_n⋯S₁.
- A `PeriodMap` struct provides `mul!`, `size` and `eltype`, applying the S's with two buffers.
- The dense case forms the product directly.

**Methods:**
- `:eigen` (dense S, n ≤ 3000): `eigen(Matrix(P))`.
- `:lu` (single sparse S): sparse LU of (S − σI), σ = 1 + √eps, then a few inverse-iteration solves from `PI.pdf.p`.
  - This is the "(S − I)p = 0 via sparse LU" idea.
  - Because S does not conserve mass exactly, λ₁ is only close to 1. A bordered normalisation system would then give a slightly different p, whereas shift-invert converges to exactly the fixed point that `advance!` iterates to.
- `:arnoldi` (large or periodic sparse): ArnoldiMethod.jl `partialschur(op; nev = 4, which = :LM)` plus `partialeigen`, using the threaded matvec.
- `:auto`: dense → `:eigen`; single sparse → `:lu`; otherwise → `:arnoldi`.

**Result handling (paper §4 caveat):**
- Pick the eigenvalue closest to 1, take the real part, and normalise with `wᵀv`, which fixes both sign and mass.
- `info = (λ, residual = ‖Sp−λp‖/‖p‖, imag_ratio, negative_mass = ∫max(−p,0), method)`.
- `@warn` if `|imag λ| > 1e-8`, `imag_ratio > 1e-6`, or `negative_mass > 1e-3`.
- Write the result to `PI.pdf.p` and set `step_idx = 0`, so a periodic PI resumes with S₁ (plus `t = 0` if `reset_t`).

**New dependency:** `ArnoldiMethod` (0.4, pure Julia, deps LinearAlgebra/Random/StaticArrays; already in the depot).

### A6. Normalisation per step via w_S = Sᵀw
- Precompute `w_S[n] = Sₙᵀ·vec(w)` once per step matrix, where `w` is the tensor product of the axis weights. Store it in a new trailing `PathIntegration` field `stepMX_wts`.
- `advance!` then computes `mass = dot(w_S[idx], vec(p))` *before* `mul!`, and scales by `1/mass`. This replaces the generic product-iterator `_integrate(p_temp, …)` (pathintegration.jl:150-164).
- Recompute `w_S` in `recompute_stepMX!`. It is `nothing` when `pre_compute = false`.

### A7. Memory for d ≥ 3
- **Relative sparsity threshold (paper Eq. 37: drop |S_ij| < max|S|·1e-8).** Add `SparseMX(; sparse_rtol = 0.0)`.
  - Rows are pre-filtered with the absolute `sparse_tol`, and each buffer tracks its max |v|.
  - If `sparse_rtol > 0`, assembly filters with τ = max(`sparse_tol`, `sparse_rtol`·max|S|) and recomputes `colptr` while copying.
  - The defaults (abs 1e-6, rel 0) stay, so current results are unchanged; the paper's setting is `SparseMX(sparse_tol = 0, sparse_rtol = 1e-8)`.
  - `sparse_rtol` is passed through `IK.kwargs` like `sparse_tol`.
- **`Int32` indices:** see A2.
- **Deferred (longer term):** applying S without storing it, i.e. recomputing rows every step. After A1 this is cheap enough to try, and since it is compute-bound it scales with threads. Not implemented here.

---

## Part B — Parallel row computation

### Per-task state (full trace through `fill_stepMX!`)
| Written per row / per quadrature node | Per task? |
|---|---|
| `IK.x1` | yes |
| `sdestep.x0`, `sdestep.x1` | yes |
| `steptracer.temp`, `.tempI` (Newton buffers) | yes |
| `method.drift.ks`, `.temp` (RK stages, aliased by `BT.a/b[].val`; **the user's `RK4()` object, shared across PIs**) | yes, via `deepcopy(method)` |
| `sdestep.t0`/`t1` Refs | own Refs; `set_t0t1!` on every workspace |
| `di.x`, `di.w`, `Q_integrate` (plus `temp`/`res` for the generic kernel); QuadGK `int_limits`, `Q_integrate`, `res`, `res0` | yes |
| `IK.temp.itpVs` plus the `idx_it`/`val_it` iterators that reference them; `itpM`; kernel workspace (`mark`/`touched`, or `Bs`/`c`/`KR`) | yes; the iterators are rebuilt |
| VIO: `SDE_VIO.ID` Ref (inside the sde), `ID_dyn`, `Q_aux`, `ti`, `xi`, `xi2`, `vitemp`, `v_i` | yes, via `deepcopy(sdestep)` |
| step matrix | dense: disjoint rows, shared; sparse: per-chunk row buffers (A2) |

Shared read-only: `sde` (f, g, `par`) for `SDE`; the Symbolics RuntimeGeneratedFunctions; `IK.pdf` (only `axes[i].xs/wts/itp` are read); `IK.t`; `IK.kwargs`. QuadGK's rule cache is lock-protected.

With A3, the result of every row depends only on the row, so every mode and thread count gives a bitwise-identical S.

### Sparse representation in parallel
Each task gets a **contiguous** chunk of rows, i.e. a column block of Sᵀ, with its own A2 row buffer. `nnz_per_row` is shared, but each i is written by exactly one task. Assembly concatenates the buffers in chunk order. No locks and no sorting are needed, and the result is exact for any chunking.

Rejected: a lock around `setindex!` (serialises writes and keeps the O(n) shifts), and COO triplets plus `sparse(I,J,V)` (3× memory and a sort).

### Types, copies, drivers
- **Types** (`src/types.jl`, exported):
  - `abstract type RowComputation`;
  - `SerialRowComputation()`;
  - `ThreadedRowComputation(N_threads = Threads.nthreads())` and `BatchRowComputation(N_threads = …)`, each with a positional or keyword `N_threads ≥ 1` (the number of chunks/workspaces).
- **`task_copy`** (new `src/integration/rowcomputation.jl`):
  - `IntegrationKernel`: via the existing outer constructor; shares `pdf`, `t`, `kwargs`.
  - `SDEStep{…,<:SDE}`: shares `sde`; copies `method`, `x0/x1`, the Refs and the tracer buffers.
  - `AbstractSDEStep` fallback: `deepcopy`.
  - `DiscreteIntegrator`: shares `x_ref/w_ref`. `QuadGKIntegrator` and `NonSmoothDiscreteIntegrator` are copied too.
  - `IK_temp`: via a new `IK_temp(itpVs, itpM, kernel)` constructor that builds the iterators; this constructor is also used at pathintegration.jl:79.
- **Row kernel:** `compute_row!(out, IK, i, idx, smart_integration)` is the current loop body. `fill_rows!` runs over a chunk using `CartesianIndices(IK.pdf.p)`. The workspace struct holds `IKs = [IK; task_copy(IK)…]`, the balanced contiguous `chunks`, and the sparse row buffers, which are reused across time intervals.
- **Drivers:**
  - Serial: one chunk on the original IK, sharing all of the above code.
  - Threaded: `Threads.@threads :static for c in eachindex(chunks)`. It indexes by chunk, not `threadid()`. When already inside a `@threads` region (`ccall(:jl_in_threaded_region, Cint, ()) != 0`), it runs the chunks sequentially, which gives the same result.
  - Batch: `Polyester.@batch per=thread for c in 1:nchunks` calling `_fill_chunk!(state, c)`. `state` is one non-bits struct, because Polyester turns captured top-level `Vector{Float64}` into `PtrArray`.
- **`fill_stepMX_ts!`:** builds the workspaces once, calls `set_t0t1!` on every workspace per interval, then `fill_stepMX!(stepMX[jₜ], ws, rc, smart_integration)`.
- **API:** `PathIntegration(…; rowcomputation = default_rowcomputation())`, with `default_rowcomputation() = Threads.nthreads() > 1 ? ThreadedRowComputation() : SerialRowComputation()`. It is stored in `IK.kwargs`, so `recompute_*` reuses it. Docstrings cover all new keywords and types.
- **Dependencies:** `Polyester` (0.7) and `ArnoldiMethod` (0.4) as hard deps in Project.toml.

---

## Implementation order
1. Baseline: `benchmark/stepmx_benchmark.jl` (committed) with S-build times for d=1/d=2 sparse and dense, and a steady-state solve. Serialise current S matrices of the test systems to the scratchpad for regression.
2. A3, then A1 + A2 (kernel/assembly refactor), A4, A6, A7. Run the tests after each.
3. Part B.
4. A5.

## Verification
- **Regression** vs the serialised baseline: sparse S within 1e-13 relative (only A3 rounding changes); dense-interpolation S within 1e-12.
- **`test/row_kernel_test.jl`:**
  - SPA vs `generic_row_kernel = true` are exactly `==`, for d=1 Cubic, d=2 Quintic, and mixed Cubic×Chebyshev, with both SparseMX and DenseMX output.
  - Dense vectorised vs generic within 1e-12 (d=1 and d=2 Chebyshev).
  - Sparse S equals dense S with entries ≤ tol zeroed, exactly.
  - `sparse_rtol` drops exactly the entries ≤ 1e-8·max.
  - `Int32` gives an S and `advance!` identical to `Int`.
- **Quadrature unit test:** rescaling a→b then c→d gives nodes bitwise equal to rescaling straight to c→d; a `NewtonCotesIntegrator` PathIntegration builds.
- **A6:** `advance!` matches `_integrate`-normalised S·p (≈1e-14); after `recompute_PI!(par=…)`, `advance!` equals a fresh PI.
- **`test/steady_state_test.jl`:**
  - `steady_state!` vs `advance_till_converged!` within 1e-5 for scalar cubic (Cubic, Chebyshev), Duffing d=2 N=21, and the periodic system;
  - `:eigen`/`:lu`/`:arnoldi` agree to 1e-8;
  - ∫p = 1, small residual, and `step_idx == 0`.
- **`test/row_computation_test.jl`:**
  - Serial / `Threaded(1,3,50)` / `Threaded()` / `Batch(3)` / `Batch()` give exactly `==` S (plus equal `rowvals`/colptr) for d=1 sparse, d=1 dense, d=2 Newton, QuadGK, and periodic with 3 intervals;
  - `recompute_PI!(par=…)` with threads equals a fresh serial PI;
  - constructor forms and `ArgumentError`;
  - `default_rowcomputation()`;
  - a nested `@threads` build doesn't error.
- **Full suite:** `julia --project -e 'using Pkg; Pkg.test()'`, and again with `JULIA_NUM_THREADS=4`.
- **Benchmark** (`julia -t 16`): serial before/after A1/A2/A4 (expect ~2× at N=81 for d=2 quintic, growing like N^d/(o+1)^d), thread scaling of Threaded/Batch for n ∈ {2,4,8,16}, and `steady_state!` vs `advance_till_converged!` time.

## Out of scope
- Dynamic load balancing (more chunks than threads).
- Parallelising over time intervals.
- Thread-safety of `PI.pdf(x)` (writes `axis.temp`).
- Several PIs sharing one `RK4()` object built concurrently.
- Applying S without storing it (A7, deferred).
- Vibro-impact: it keeps the generic kernel and is still broken under Symbolics 7.
