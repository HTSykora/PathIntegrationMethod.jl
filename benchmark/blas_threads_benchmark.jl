# Effect of single_threaded_blas (a single BLAS thread during the parallel step matrix computation) with dense interpolations
# The configurations are timed interleaved, 40 times each (the timings of ~10 ms are noisy)
# Run from the package folder with: julia --project -t auto benchmark/blas_threads_benchmark.jl
using PathIntegrationMethod, Statistics
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]
for N in (41, 61)
    PI = PathIntegration(SDE((f1, f2), g2, [0.5,0.25,1.0]), RK4(), 0.01, ChebyshevAxis(-4.,4.,N), ChebyshevAxis(-4.,4.,N); rowcomputation = SerialRowComputation())
    configs = Tuple{String,Any,Bool}[("serial", SerialRowComputation(), false)]
    for n in unique([Threads.nthreads() ÷ 2, Threads.nthreads()]), (name, RC) in (("Threaded", ThreadedRowComputation), ("Batch", BatchRowComputation)), blas1 in (false, true)
        push!(configs, ("$name($n)$(blas1 ? ", single_threaded_blas" : "")", RC(n; single_threaded_blas = blas1), blas1))
    end
    times = Dict(lbl => Float64[] for (lbl, _, _) in configs)
    for (lbl, rc, _) in configs; recompute_stepMX!(PI; rowcomputation = rc); end
    for _ in 1:40, (lbl, rc, _) in configs # interleaved
        push!(times[lbl], @elapsed recompute_stepMX!(PI; rowcomputation = rc))
    end
    println("d=2 Chebyshev N=$N ($(length(PI.pdf)) rows), $(BLAS.get_num_threads()) BLAS threads; recompute_stepMX! median (min) of 40:")
    for (lbl, _, _) in configs
        println("   ", rpad(lbl, 36), round(1e3median(times[lbl]), sigdigits = 3), " ms (", round(1e3minimum(times[lbl]), sigdigits = 3), " ms)")
    end
end
