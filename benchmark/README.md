# Benchmarks

Scripts for measuring the performance of the step matrix computation. Run them from the package folder.

- `stepmx_benchmark.jl`: serial timings for d = 1 and d = 2 systems with sparse and dense interpolations: step matrix computation (as in `PathIntegration`), `recompute_stepMX!`, number of nonzeros, one `advance!`, `advance_till_converged!` and `steady_state!`.
  Run with `julia --project -t 1 benchmark/stepmx_benchmark.jl`.
- `parallel_benchmark.jl`: thread scaling of the step matrix computation with `ThreadedRowComputation(n)` and `BatchRowComputation(n)` compared with `SerialRowComputation()` (it also checks that the step matrices are identical).
  Run with `julia --project -t auto benchmark/parallel_benchmark.jl`.
- `regression_snapshot.jl`: stores the step matrices of a set of test systems (`save`) or compares the current ones with the stored ones (`compare`), to check that a performance change does not change the results.
  Run with `julia --project benchmark/regression_snapshot.jl save|compare`.

`results/` holds the outputs (AMD Ryzen 7 PRO 7840U, 8 cores / 16 threads):

- `stepmx_benchmark_before.txt`: commit d5c342b (before the optimisations)
- `stepmx_benchmark_after_partA.txt`: commit 9f51dfa (sparse accumulator, direct CSC assembly, stateless quadrature, vectorised dense rows, `steady_state!`). Its `recompute_stepMX!` times are still dominated by resetting the old step matrix element by element.
- `stepmx_benchmark_after.txt`: after the parallelisation changes (type-stable Runge–Kutta stages, fast reset of the step matrix), serial
- `parallel_benchmark.txt`: thread scaling (with 8 BLAS threads; with dense interpolations `BLAS.set_num_threads(1)` is about 2× faster at 16 threads)
- `regression_after_partA.txt`: comparison with the step matrices of d5c342b

The `.jls` snapshot files are not committed.
