# Timing of the step matrix computation and of the steady-state solution
# Run from the package folder with: julia --project -t auto benchmark/stepmx_benchmark.jl
using PathIntegrationMethod
const PIM = PathIntegrationMethod

# Scalar system: dx = (x - x³) dt + √2 dW
f(x,p,t) = x[1] - x[1]^3
g(x,p,t) = sqrt(2)
# Duffing oscillator: dx = v dt, dv = (-2ζv + x - λx³) dt + σ dW
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]
par = [0.5, 0.25, 1.0] # ζ, λ, σ

scalar(axis, N; kwargs...) = PathIntegration(SDE(f, g), RK4(), 0.01, axis(-3., 3., N); kwargs...)
duffing(axis, N; kwargs...) = PathIntegration(SDE((f1, f2), g2, copy(par)), RK4(), 0.01, axis(-4., 4., N), axis(-4., 4., N); kwargs...)

# (label, constructor, N)
cases = (("d=1 cubic", N -> scalar(CubicAxis, N), 101),
         ("d=1 Chebyshev", N -> scalar(ChebyshevAxis, N), 61),
         ("d=2 quintic", N -> duffing(QuinticAxis, N), 41),
         ("d=2 quintic", N -> duffing(QuinticAxis, N), 81),
         ("d=2 Chebyshev", N -> duffing(ChebyshevAxis, N), 31))

# Minimum time of `n` evaluations of `f()`
besttime(f, n = 3) = minimum(@elapsed(f()) for _ in 1:n)

println("Julia threads: ", Threads.nthreads())
for (lbl, PI_N, N) in cases
    PI_small = PI_N(21) # compilation
    advance_till_converged!(PI_small; Tmax = 1.0)
    # (the coarse warm-up grids trigger accuracy warnings)
    isdefined(PIM, :steady_state!) && Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        foreach(m -> PIM.steady_state!(PI_small; method = m), (:auto, :lu))
    end
    PI = PI_N(N)
    t_S = besttime(() -> recompute_stepMX!(PI))
    S = PI.stepMX[1]
    nz = S isa AbstractMatrix{<:Number} && !(S isa Matrix) ? length(PIM.SparseArrays.nonzeros(parent(S))) : length(S)
    t_step = besttime(() -> advance!(PI), 100)
    t_conv = besttime(() -> (reinit_PI_pdf!(PI); advance_till_converged!(PI; Tmax = 20.0)))
    line = "$lbl N = $N: S build $(round(t_S, sigdigits = 3)) s, nnz(S) = $nz, one step $(round(1e6t_step, sigdigits = 3)) μs, advance_till_converged! $(round(t_conv, sigdigits = 3)) s"
    if isdefined(PIM, :steady_state!)
        t_ss = besttime(() -> (reinit_PI_pdf!(PI); PIM.steady_state!(PI)))
        line *= ", steady_state! $(round(t_ss, sigdigits = 3)) s"
    end
    println(line)
end
