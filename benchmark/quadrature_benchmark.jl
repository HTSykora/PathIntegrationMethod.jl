# Quadrature rules of the Chapman–Kolmogorov integral (GaussLegendreIntegrator, GaussHermiteIntegrator) with noise on one and on several coordinates:
# the difference of the stationary PDF to a reference quadrature with many nodes, and for noise on several coordinates the error to the exact stationary PDF.
# Run from the package folder with: julia --project -t auto benchmark/quadrature_benchmark.jl
using PathIntegrationMethod, LinearAlgebra, Printf
const PIM = PathIntegrationMethod
quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
stationary(PI) = (quiet(() -> steady_state!(PI)); PI)

println("Julia threads: ", Threads.nthreads())
# Noise on the last coordinate: Duffing oscillator dx = v dt, dv = (-2ζv + x - λx³) dt + σ dW
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]
duffing(Δt, di, bt) = PathIntegration(SDE((f1, f2), g2, [0.5, 0.25, 1.0]), RK4(), Δt, QuinticAxis(-4., 4., 41), QuinticAxis(-4., 4., 41); discreteintegrator = di, backtracing = bt)
Ns = (3, 5, 7, 9, 11, 15, 21, 31)
println("\nDuffing oscillator (quintic 41², RK4): L¹ difference of the stationary PDF to Gauss–Legendre with 121 nodes")
for bt in (NewtonBacktracing(), ExplicitBacktracing()), Δt in (0.01, 0.05, 0.2)
    ref = stationary(duffing(Δt, GaussLegendreIntegrator(121), bt))
    for (lbl, rule) in (("Gauss–Legendre", GaussLegendreIntegrator), ("Gauss–Hermite", GaussHermiteIntegrator))
        @printf("  %-22s Δt = %-5s %-15s", bt, Δt, lbl)
        foreach(N -> @printf(" %2d: %.0e", N, integrate_diff(stationary(duffing(Δt, rule(N), bt)).pdf, ref.pdf.p)), Ns)
        println()
    end
end

# Noise on two coordinates: dx = v dt, dv = (-x - 2ζv + y) dt + σ₁ dW₁, dy = -αy dt + σ₂ dW₂, stationary PDF N(0, P), A P + P Aᵀ + B Bᵀ = 0
ζ, α, σ₁, σ₂ = 0.5, 1.0, 0.5, 0.7
h1(x,p,t) = x[2]
h2(x,p,t) = -x[1] - 2p[1]*x[2] + x[3]
h3(x,p,t) = -p[2]*x[3]
s2(x,p,t) = p[3]
s3(x,p,t) = p[4]
P = lyap([0 1 0; -1 -2ζ 1; 0 0 -α], Diagonal([0, σ₁^2, σ₂^2]))
axes3(N) = Tuple(QuinticAxis(-L, L, N) for L in 4.5 .* sqrt.(diag(P)))
oscillator(Δt, N; kwargs...) = PathIntegration(SDE((h1, h2, h3), (s2, s3), [ζ, α, σ₁, σ₂]), RK4(), Δt, axes3(N)...; backtracing = ExplicitBacktracing(), kwargs...)
function exact_error(PI)
    PIe = PathIntegration(PI.step_dynamics, PI.ts, PI.pdf.axes...; f_init = (x...) -> exp(-0.5*dot(collect(x), P \ collect(x))), pre_compute = false)
    PIe.pdf.p ./= integrate(PIe.pdf)
    integrate_diff(PI.pdf, PIe.pdf)
end
println("\nOscillator driven by white noise and an Ornstein–Uhlenbeck process (noise on 2 coordinates, quintic 25³, RK4, ExplicitBacktracing)")
println("  L¹ error to the exact stationary PDF (default Gauss–Hermite 7²):")
for Δt in (0.1, 0.05, 0.025)
    @printf("    Δt = %-6s %.2e\n", Δt, exact_error(stationary(oscillator(Δt, 25))))
end
println("  Δt = 0.05: L¹ difference of the stationary PDF to Gauss–Legendre 41², and the time of the step matrix computation:")
ref = stationary(oscillator(0.05, 25; discreteintegrator = GaussLegendreIntegrator(41, dim = 2)))
for (lbl, di) in (("Gauss–Hermite 3²", GaussHermiteIntegrator(3, dim = 2)), ("Gauss–Hermite 5²", GaussHermiteIntegrator(5, dim = 2)),
                  ("Gauss–Hermite 7²", GaussHermiteIntegrator(7, dim = 2)), ("Gauss–Hermite 11²", GaussHermiteIntegrator(11, dim = 2)),
                  ("Gauss–Legendre 11²", GaussLegendreIntegrator(11, dim = 2)), ("Gauss–Legendre 15²", GaussLegendreIntegrator(15, dim = 2)),
                  ("Gauss–Legendre 21²", GaussLegendreIntegrator(21, dim = 2)), ("Gauss–Legendre 31²", GaussLegendreIntegrator(31, dim = 2)))
    PI = oscillator(0.05, 25; discreteintegrator = di)
    t = minimum(@elapsed(recompute_stepMX!(PI)) for _ in 1:2)
    @printf("    %-20s %4d nodes: %.1e   step matrix %.2f s\n", lbl, length(PI.IK.discreteintegrator.x), integrate_diff(stationary(PI).pdf, ref.pdf.p), t)
end
