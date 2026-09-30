# Benchmarks

Scripts for measuring the performance of the step matrix computation. Run them from the package folder.

- `stepmx_benchmark.jl`: serial timings for d = 1 and d = 2 systems with sparse and dense interpolations: step matrix computation (as in `PathIntegration`), `recompute_stepMX!`, number of nonzeros, one `advance!`, `advance_till_converged!` and `steady_state!`.
  Run with `julia --project -t 1 benchmark/stepmx_benchmark.jl`.
- `parallel_benchmark.jl`: thread scaling of the step matrix computation with `ThreadedRowComputation(n)` and `BatchRowComputation(n)` compared with `SerialRowComputation()` (it also checks that the step matrices are identical).
  Run with `julia --project -t auto benchmark/parallel_benchmark.jl`.
- `dense_layout_benchmark.jl`: dense step matrices stored as `S` or as `transpose(Sᵀ)`: matrix–vector products (BLAS gemv 'N' vs 'T') and end-to-end `advance!`, `advance_till_converged!` and `steady_state!` timings.
  Run with `julia --project -t auto benchmark/dense_layout_benchmark.jl`.
- `dense_kernel_benchmark.jl`: dense interpolations, the tensor row kernel compared with the old generic row kernel, with 1 … `Threads.nthreads()` threads.
  Run with `julia --project -t auto benchmark/dense_kernel_benchmark.jl`.
- `blas_threads_benchmark.jl`: effect of `single_threaded_blas = true` in `ThreadedRowComputation`/`BatchRowComputation` with dense interpolations (interleaved repetitions, medians).
  Run with `julia --project -t auto benchmark/blas_threads_benchmark.jl`.
- `regression_snapshot.jl`: stores the step matrices of a set of test systems (`save`) or compares the current ones with the stored ones (`compare`), to check that a performance change does not change the results.
  Run with `julia --project benchmark/regression_snapshot.jl save|compare`.
- `duffing_2d_profile.jl`: Duffing oscillator (d = 2) with uniform (quintic, Chebyshev) and mixed dense-sparse axes: timings and allocations of the step matrix computation, `advance!`, `integrate`, `integrate_diff`, the evaluation `PI(x, v)` and `update_mPDFs!` (`bench`), and a snapshot of the results (step matrices, PDFs after 20 steps, integrals, values, marginal PDFs) to check that a change does not change them (`save`, `compare`).
  Run with `julia --project -t auto benchmark/duffing_2d_profile.jl bench|save|compare`.
- `duffing_2d_ab.jl`: the same Duffing oscillator, also with mixed sparse-sparse and dense-dense axes, to compare two versions of the package (run alternately with each version, see the script).
  Run with `julia --project -t auto benchmark/duffing_2d_ab.jl`.
- `backtracing_benchmark.jl`: `NewtonBacktracing`, `ExplicitBacktracing` and `StrangSplitting` on the Duffing oscillator: the time to the first step matrix of a new SDE, the step matrix computation (serial and threaded), and the error of the stationary PDF (exact solution) for several time steps.
  Run with `julia --project -t auto benchmark/backtracing_benchmark.jl`.

`results/` holds the outputs (AMD Ryzen 7 PRO 7840U, 8 cores / 16 threads):

- `stepmx_benchmark_before.txt`: commit d5c342b (before the optimisations)
- `stepmx_benchmark_after_partA.txt`: commit 9f51dfa (sparse accumulator, direct CSC assembly, stateless quadrature, vectorised dense rows, `steady_state!`). Its `recompute_stepMX!` times are still dominated by resetting the old step matrix element by element.
- `stepmx_benchmark_after.txt`: after the parallelisation changes (type-stable Runge–Kutta stages, fast reset of the step matrix, dense step matrices stored as `transpose(Sᵀ)`), serial
- `parallel_benchmark.txt`: thread scaling, also with `single_threaded_blas = true` for dense interpolations
- `dense_layout_benchmark.txt`: after storing the dense step matrices as `transpose(Sᵀ)`; `dense_layout_before.txt`: the same end-to-end timings with `S` stored directly
- `dense_kernel_benchmark.txt`, `blas_threads_benchmark.txt`
- `regression_after_partA.txt`: comparison with the step matrices of d5c342b
- `duffing_2d_profile_before.txt`, `duffing_2d_profile_after.txt`: `duffing_2d_profile.jl bench` before (commit 381602e) and after the type stability fixes (loops over axes of different types, `integrate`, the evaluation `PI(x, v)`)
- `duffing_2d_ab.txt`: `duffing_2d_ab.jl` with commit 381602e (old) and after the type stability fixes (new), two alternating runs each
- `backtracing_benchmark.txt`: `backtracing_benchmark.jl`

The `.jls` snapshot files are not committed.
