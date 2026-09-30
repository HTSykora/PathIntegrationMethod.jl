# Back-tracing methods (NewtonBacktracing, ExplicitBacktracing, StrangSplitting) on the Duffing oscillator
# dx = v dt, dv = (-2ζv + x - λx³) dt + σ dW, which has the exact stationary PDF ∝ exp(-(4ζ/σ²)(v²/2 - x²/2 + λx⁴/4)):
# time to the first step matrix of a new SDE, step matrix computation, and the error of the stationary PDF.
# Run from the package folder with: julia --project -t auto benchmark/backtracing_benchmark.jl
using PathIntegrationMethod, Printf, Logging
const PIM = PathIntegrationMethod

const ζ, λ, σ = 0.5, 0.25, 1.0
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]
p_exact(x, v) = exp(-(4ζ/σ^2)*(v^2/2 - x^2/2 + λ*x^4/4))
duffing() = SDE((f1, f2), g2, [ζ, λ, σ])
bts = (NewtonBacktracing(), ExplicitBacktracing(), StrangSplitting())
besttime(f, n = 5) = minimum(@elapsed(f()) for _ in 1:n)

println("Julia threads: ", Threads.nthreads())
# a new drift (a new function) needs a new symbolic Newton step (NewtonBacktracing), or the compilation of the explicit steps
for bt in bts
    PathIntegration(duffing(), RK4(), 0.05, QuinticAxis(-4., 4., 21), QuinticAxis(-4., 4., 21); backtracing = bt) # compilation of the package code
end
new_drifts = [@eval((x,p,t) -> -2p[1]*x[2] + x[1] - p[2]*x[1]^3 + $(1e-3i)*x[1]) for i in 1:3length(bts)]
for (j, bt) in enumerate(bts)
    ts = map(1:3) do i
        h2 = new_drifts[3(j - 1) + i]
        @elapsed PathIntegration(SDE((f1, h2), g2, [ζ, λ, σ]), RK4(), 0.05, QuinticAxis(-4., 4., 41), QuinticAxis(-4., 4., 41); backtracing = bt)
    end
    @printf("%-22s first step matrix of a new SDE (quintic 41 × 41, RK4): %.2f s\n", bt, minimum(ts))
end

println("\nStep matrix computation (RK4, Δt = 0.05)")
for (lbl, axes) in (("quintic 41 × 41", () -> (QuinticAxis(-4., 4., 41), QuinticAxis(-4., 4., 41))),
                    ("quintic 81 × 81", () -> (QuinticAxis(-4., 4., 81), QuinticAxis(-4., 4., 81))),
                    ("Chebyshev 25 × 25", () -> (ChebyshevAxis(-4., 4., 25), ChebyshevAxis(-4., 4., 25))))
    for bt in bts
        PI = PathIntegration(duffing(), RK4(), 0.05, axes()...; backtracing = bt)
        t_s = besttime(() -> recompute_stepMX!(PI; rowcomputation = SerialRowComputation()))
        t_p = besttime(() -> recompute_stepMX!(PI; rowcomputation = ThreadedRowComputation()))
        S = PIM.storage_matrix(PI.stepMX[1])
        t_a = besttime(() -> advance!(PI), 100)
        @printf("  %-18s %-22s serial %7.1f ms | %2d threads %6.1f ms | nnz(S) %7d | advance! %6.1f μs\n", lbl, bt, 1e3t_s, Threads.nthreads(), 1e3t_p,
            S isa PIM.SparseMatrixCSC ? PIM.nnz(S) : length(S), 1e6t_a)
    end
end

function stationary_error(bt, Δt, N)
    axes() = (QuinticAxis(-4., 4., N), QuinticAxis(-4., 4., N))
    PI = PathIntegration(duffing(), RK4(), Δt, axes()...; backtracing = bt)
    with_logger(NullLogger()) do
        steady_state!(PI)
    end
    PIe = PathIntegration(duffing(), RK4(), Δt, axes()...; f_init = p_exact, pre_compute = false)
    PIe.pdf.p ./= integrate(PIe.pdf)
    integrate_diff(PI.pdf, PIe.pdf)
end
Δts = (0.1, 0.05, 0.025, 0.0125)
for (N, grid_bts) in ((81, bts), (121, (StrangSplitting(),)))
    println("\nL¹ error of the stationary PDF (RK4, quintic $N × $N)")
    @printf("  %-22s", "Δt"); foreach(Δt -> @printf("%10s", Δt), Δts); println("   orders")
    for bt in grid_bts
        errs = [stationary_error(bt, Δt, N) for Δt in Δts]
        @printf("  %-22s", bt); foreach(e -> @printf("%10.2e", e), errs)
        println("   ", join((@sprintf("%.2f", log2(errs[i] / errs[i+1])) for i in 1:length(errs)-1), ", "))
    end
end
println("(at small Δt the error of StrangSplitting is limited by the interpolation error of the grid, which grows like h⁶/Δt)")
