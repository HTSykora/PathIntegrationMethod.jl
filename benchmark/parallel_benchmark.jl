# Thread scaling of the step matrix computation (SerialRowComputation vs ThreadedRowComputation and BatchRowComputation)
# Run from the package folder with: julia --project -t auto benchmark/parallel_benchmark.jl
using PathIntegrationMethod

# Duffing oscillator: dx = v dt, dv = (-2ζv + x - λx³) dt + σ dW
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]
par = [0.5, 0.25, 1.0] # ζ, λ, σ
# Linear system with d = 3
h1(x,p,t) = x[2]
h2(x,p,t) = -x[1] - x[2] + x[3]
h3(x,p,t) = -x[3] - x[2]
g3(x,p,t) = 0.8

duffing(axis, N) = PathIntegration(SDE((f1, f2), g2, copy(par)), RK4(), 0.01, axis(-4., 4., N), axis(-4., 4., N); rowcomputation = SerialRowComputation())
linear3(axis, N) = PathIntegration(SDE((h1, h2, h3), g3), RK4(), 0.05, axis(-3., 3., N), axis(-3., 3., N), axis(-3., 3., N); rowcomputation = SerialRowComputation())

# (label, constructor, N)
cases = (("d=2 quintic", N -> duffing(QuinticAxis, N), 81),
         ("d=2 Chebyshev", N -> duffing(ChebyshevAxis, N), 41),
         ("d=3 cubic", N -> linear3(CubicAxis, N), 21))

# Minimum time of `n` evaluations of `f()`
besttime(f, n = 3) = minimum(@elapsed(f()) for _ in 1:n)
S_matrices(PI) = [Matrix(S) for S in PI.stepMX]

nthreads = Threads.nthreads()
ns = unique([filter(≤(nthreads), [1, 2, 4, 8, 16, 32]); nthreads])
println("Julia threads: ", nthreads, ", BLAS threads: ", PathIntegrationMethod.BLAS.get_num_threads())
for (lbl, PI_N, N) in cases
    for rc in (SerialRowComputation(), ThreadedRowComputation(2), BatchRowComputation(2)) # compilation
        recompute_stepMX!(PI_N(7); rowcomputation = rc)
    end
    PI = PI_N(N)
    t_serial = besttime(() -> recompute_stepMX!(PI; rowcomputation = SerialRowComputation()))
    S_serial = S_matrices(PI)
    println("$lbl N = $N ($(length(PI.pdf)) rows): serial $(round(t_serial, sigdigits = 3)) s")
    for (name, RC) in (("Threaded", ThreadedRowComputation), ("Batch", BatchRowComputation))
        line = "    $(rpad(name, 8))"
        for n in ns
            t = besttime(() -> recompute_stepMX!(PI; rowcomputation = RC(n)))
            S_matrices(PI) == S_serial || error("$name($n) gives a different step matrix")
            line *= "  n = $n: $(round(t, sigdigits = 3)) s ($(round(t_serial / t, digits = 1))×)"
        end
        println(line)
    end
end
