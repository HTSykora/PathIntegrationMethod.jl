# Benchmarks

Scripts for measuring the performance of the step matrix computation. Run them from the package folder.

- `stepmx_benchmark.jl`: S build time, number of nonzeros, time of one `advance!`, `advance_till_converged!` and `steady_state!` for d = 1 and d = 2 systems with sparse and dense interpolations.
  Run with `julia --project -t 1 benchmark/stepmx_benchmark.jl` (serial) or `-t auto`.
- `regression_snapshot.jl`: stores the step matrices of a set of test systems (`save`) or compares the current ones with the stored ones (`compare`), to check that a performance change does not change the results.
  Run with `julia --project benchmark/regression_snapshot.jl save|compare`.

`results/` holds the outputs: `stepmx_benchmark_before.txt` (commit d5c342b) and `stepmx_benchmark_after_partA.txt` (sparse accumulator, direct CSC assembly, stateless quadrature, vectorised dense rows, `steady_state!`), and `regression_after_partA.txt` (comparison with the d5c342b step matrices). The `.jls` snapshot files are not committed.
