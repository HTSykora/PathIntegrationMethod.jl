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

`results/` holds the outputs (AMD Ryzen 7 PRO 7840U, 8 cores / 16 threads):

- `stepmx_benchmark_before.txt`: commit d5c342b (before the optimisations)
- `stepmx_benchmark_after_partA.txt`: commit 9f51dfa (sparse accumulator, direct CSC assembly, stateless quadrature, vectorised dense rows, `steady_state!`). Its `recompute_stepMX!` times are still dominated by resetting the old step matrix element by element.
- `stepmx_benchmark_after.txt`: after the parallelisation changes (type-stable Runge–Kutta stages, fast reset of the step matrix, dense step matrices stored as `transpose(Sᵀ)`), serial
- `parallel_benchmark.txt`: thread scaling, also with `single_threaded_blas = true` for dense interpolations
- `dense_layout_benchmark.txt`: after storing the dense step matrices as `transpose(Sᵀ)`; `dense_layout_before.txt`: the same end-to-end timings with `S` stored directly
- `dense_kernel_benchmark.txt`, `blas_threads_benchmark.txt`
- `regression_after_partA.txt`: comparison with the step matrices of d5c342b

The `.jls` snapshot files are not committed.
