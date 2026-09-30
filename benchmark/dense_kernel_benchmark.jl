# Dense interpolations: the tensor row kernel (one BLAS matrix product per row) compared with the old generic row kernel
# (the DiscreteIntegrator summing the full row at every quadrature node), with 1 … Threads.nthreads() threads
# Run from the package folder with: julia --project -t auto benchmark/dense_kernel_benchmark.jl
using PathIntegrationMethod
f(x,p,t) = x[1] - x[1]^3
g(x,p,t) = sqrt(2)
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]
duffing(N; kw...) = PathIntegration(SDE((f1, f2), g2, [0.5,0.25,1.0]), RK4(), 0.01, ChebyshevAxis(-4.,4.,N), ChebyshevAxis(-4.,4.,N); rowcomputation = SerialRowComputation(), kw...)
scalar(N; kw...) = PathIntegration(SDE(f, g), RK4(), 0.01, ChebyshevAxis(-3.,3.,N); rowcomputation = SerialRowComputation(), kw...)
besttime(f, n = 5) = minimum(@elapsed(f()) for _ in 1:n)
ns = filter(≤(Threads.nthreads()), [1, 2, 4, 8, 16, 32])
println("Julia threads: ", Threads.nthreads(), ", BLAS threads: ", BLAS.get_num_threads())
for (lbl, mk) in (("d=2 Chebyshev N=41", (; kw...) -> duffing(41; kw...)), ("d=2 Chebyshev N=61", (; kw...) -> duffing(61; kw...)), ("d=1 Chebyshev N=101", (; kw...) -> scalar(101; kw...)))
    PI_t = mk(); PI_g = mk(; generic_row_kernel = true) # (generic_row_kernel: internal keyword for tests and benchmarks)
    S_ref = Matrix(PI_t.stepMX[1])
    println(lbl, " (", length(PI_t.pdf), " rows)")
    for (variant, PI) in (("tensor kernel (current)", PI_t), ("old generic kernel", PI_g))
        recompute_stepMX!(PI; rowcomputation = ThreadedRowComputation(2)) # compilation
        ts = [besttime(() -> recompute_stepMX!(PI; rowcomputation = n == 0 ? SerialRowComputation() : ThreadedRowComputation(n))) for n in [0; ns]]
        err = maximum(abs, Matrix(PI.stepMX[1]) - S_ref) / maximum(abs, S_ref)
        println("  ", rpad(variant, 25), "serial ", rpad(round(1000ts[1], sigdigits = 3), 6), " ms | ",
                join(("$n: $(round(1000t, sigdigits = 3)) ms" for (n, t) in zip(ns, ts[2:end])), ", "), "   (relative difference $(round(err, sigdigits = 2)))")
    end
end
