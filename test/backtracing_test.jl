module BacktracingTest
using PathIntegrationMethod, Test
const PIM = PathIntegrationMethod

# Scalar system: dx = (x - x³) dt + √2 dW, stationary PDF ∝ exp(x²/2 - x⁴/4)
f(x,p,t) = x[1] - x[1]^3
g(x,p,t) = sqrt(2)
p_scalar(x) = exp(x^2/2 - x^4/4)
gm(x,p,t) = 0.5 + 0.5x[1]^2 # multiplicative noise
# Duffing oscillator: dx = v dt, dv = (-2ζv + x - λx³) dt + σ dW, stationary PDF ∝ exp(-(4ζ/σ²)(v²/2 - x²/2 + λx⁴/4))
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]
duffing(par = [0.5, 0.25, 1.0]) = SDE((f1, f2), g2, par)
p_duffing(x, v) = exp(-2*(v^2/2 - x^2/2 + 0.25*x^4/4))

# L¹ error of the stationary PDF
function stationary_error(PI, p_exact)
    Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        steady_state!(PI)
    end
    PI_exact = PathIntegration(PI.step_dynamics, PI.ts, PI.pdf.axes...; f_init = p_exact, pre_compute = false)
    PI_exact.pdf.p ./= integrate(PI_exact.pdf)
    integrate_diff(PI.pdf, PI_exact.pdf)
end
stepmatrices(PI) = [Matrix(S) for S in PI.stepMX]

@testset "ExplicitBacktracing gives the results of NewtonBacktracing (RK4)" begin
    for (new_PI, p_exact) in ((bt -> PathIntegration(SDE(f, g), RK4(), 0.025, ChebyshevAxis(-3., 3., 41); backtracing = bt), p_scalar),
                              (bt -> PathIntegration(duffing(), RK4(), 0.05, QuinticAxis(-4., 4., 31), QuinticAxis(-4., 4., 31); backtracing = bt), p_duffing))
        e_newton = stationary_error(new_PI(NewtonBacktracing()), p_exact)
        e_explicit = stationary_error(new_PI(ExplicitBacktracing()), p_exact)
        @test abs(e_explicit - e_newton) < 0.02e_newton
    end
end

@testset "StrangSplitting is second order" begin
    new_PI(bt, Δt) = PathIntegration(SDE(f, g), RK4(), Δt, ChebyshevAxis(-3., 3., 41); backtracing = bt)
    e1, e2 = stationary_error(new_PI(StrangSplitting(), 0.025), p_scalar), stationary_error(new_PI(StrangSplitting(), 0.0125), p_scalar)
    @test e1 / e2 > 3.5
    @test e2 < stationary_error(new_PI(NewtonBacktracing(), 0.0125), p_scalar) / 20
    @test stationary_error(PathIntegration(duffing(), RK4(), 0.05, QuinticAxis(-4., 4., 41), QuinticAxis(-4., 4., 41); backtracing = StrangSplitting()), p_duffing) < 0.003
end

# dx = (p₁ + p₂t + p₃t² + p₄t³) dt + σ dW: one step shifts the mean by ∫f(t)dt (exact up to degree 0, 1, 3 for Euler, RK2, RK4)
ft(x,p,t) = p[1] + p[2]*t + p[3]*t^2 + p[4]*t^3
F(t, p) = p[1]*t + p[2]*t^2/2 + p[3]*t^3/3 + p[4]*t^4/4
@testset "Time-dependent drift: $bt, $lbl" for bt in (ExplicitBacktracing(), StrangSplitting()),
        (lbl, method, par) in (("Euler", Euler(), [0.5, 0.0, 0.0, 0.0]), ("RK2", RK2(), [0.5, 2.0, 0.0, 0.0]), ("RK4", RK4(), [0.5, 2.0, -1.0, 3.0]))
    t0, t1, x0 = 1.0, 1.1, 0.3
    # (a fine grid: the interpolation error of the initial PDF shifts the mean by up to 6e-7 with 61 nodes)
    PI = PathIntegration(SDE(ft, (x,p,t) -> 0.5, par), method, [t0, t1], ChebyshevAxis(-2., 4., 121); f_init = x -> exp(-(x - x0)^2 / (2*0.2^2)), backtracing = bt)
    advance!(PI)
    ax = PI.pdf.axes[1]
    @test sum(ax.wts .* ax.xs .* PI.pdf.p) / sum(ax.wts .* PI.pdf.p) ≈ x0 + F(t1, par) - F(t0, par) atol = 1e-8
end

@testset "Time-periodic system: $bt" for bt in (ExplicitBacktracing(), StrangSplitting())
    fp(x,p,t) = x[1] - x[1]^3 + p[1]*sin(2π*t/p[2])
    new_PI(bt) = PathIntegration(SDE(fp, g, [1.0, 0.5]), RK4(), collect(range(0, 0.5, length = 11)), QuinticAxis(-3., 3., 61); backtracing = bt)
    PI, PI_newton = new_PI(bt), new_PI(NewtonBacktracing())
    @test length(PI.stepMX) == 10
    for _ in 1:10
        advance!(PI); advance!(PI_newton)
    end
    @test integrate_diff(PI.pdf, PI_newton.pdf.p) < (bt isa ExplicitBacktracing ? 0.01 : 0.05) # (StrangSplitting: the error of NewtonBacktracing)
end

@testset "Multiplicative noise" begin
    new_PI(bt) = PathIntegration(SDE(f, gm), RK4(), 0.02, QuinticAxis(-3., 3., 61); backtracing = bt)
    PI, PI_newton = new_PI(ExplicitBacktracing()), new_PI(NewtonBacktracing())
    steady_state!(PI); steady_state!(PI_newton)
    @test integrate_diff(PI.pdf, PI_newton.pdf.p) < 0.01
end

@testset "Drift that branches on the state (not traceable by Symbolics): $bt" for bt in (ExplicitBacktracing(), StrangSplitting())
    fb(x,p,t) = x[1] > 0 ? -x[1] : -2x[1]
    PI = PathIntegration(SDE(fb, g), RK4(), 0.05, QuinticAxis(-5., 5., 101); backtracing = bt)
    _, info = steady_state!(PI)
    @test info.residual < 1e-10
    # stationary PDF ∝ exp(-x²/2) for x > 0 and exp(-x²) for x < 0
    PI_exact = PathIntegration(PI.step_dynamics, PI.ts, PI.pdf.axes...; f_init = x -> exp(-(x > 0 ? 1 : 2)*x^2/2), pre_compute = false)
    PI_exact.pdf.p ./= integrate(PI_exact.pdf)
    @test integrate_diff(PI.pdf, PI_exact.pdf) < 0.05
end

@testset "Step matrix representations and interpolations: $bt" for bt in (ExplicitBacktracing(), StrangSplitting())
    new_PI(axes; kwargs...) = PathIntegration(duffing(), RK4(), 0.05, axes...; backtracing = bt, kwargs...)
    quintic = (QuinticAxis(-4., 4., 21), QuinticAxis(-4., 4., 21))
    S = stepmatrices(new_PI(quintic))[1]
    @test maximum(abs, stepmatrices(new_PI(quintic; stepMXtype = DenseMX()))[1] - S) ≤ 1e-6
    @test PIM.storage_matrix(new_PI(quintic; stepMXtype = DenseMX()).stepMX[1]) isa Matrix
    PI32 = new_PI(quintic; stepMXtype = SparseMX(index_type = Int32))
    @test PIM.storage_matrix(PI32.stepMX[1]) isa PIM.SparseMatrixCSC{Float64,Int32}
    @test stepmatrices(PI32)[1] == S
    S_rtol = stepmatrices(new_PI(quintic; stepMXtype = SparseMX(sparse_tol = 0.0, sparse_rtol = 1e-3)))[1]
    @test all(v -> v == 0 || abs(v) > 1e-3 * maximum(abs, S_rtol), S_rtol)
    # dense, mixed and QuadGK
    for (axes, kwargs) in (((ChebyshevAxis(-4., 4., 15), ChebyshevAxis(-4., 4., 15)), (;)),
                           ((QuinticAxis(-4., 4., 21), ChebyshevAxis(-4., 4., 15)), (;)),
                           ((CubicAxis(-4., 4., 15), CubicAxis(-4., 4., 15)), (; discreteintegrator = QuadGKIntegrator())))
        PI = new_PI(axes; kwargs...)
        advance!(PI)
        @test integrate(PI.pdf) ≈ 1
        @test all(isfinite, PI.pdf.p)
    end
end

@testset "Row computations and recompute_PI!: $bt" for bt in (ExplicitBacktracing(), StrangSplitting())
    new_PI(par, rc) = PathIntegration(duffing(par), RK4(), 0.05, QuinticAxis(-4., 4., 21), ChebyshevAxis(-4., 4., 15); backtracing = bt, rowcomputation = rc)
    par = [0.5, 0.25, 1.0]
    S = stepmatrices(new_PI(copy(par), SerialRowComputation()))
    @test stepmatrices(new_PI(copy(par), ThreadedRowComputation(3))) == S
    @test stepmatrices(new_PI(copy(par), BatchRowComputation(3))) == S
    PI = new_PI(copy(par), ThreadedRowComputation(3))
    recompute_PI!(PI; par = 0.8 .* par)
    @test stepmatrices(PI) == stepmatrices(new_PI(0.8 .* par, SerialRowComputation()))
    recompute_stepMX!(PI; t = [0.0, 0.05, 0.1]) # two time intervals
    @test length(PI.stepMX) == 2
    @test stepmatrices(PI)[1] ≈ stepmatrices(PI)[2]
end

# (measured inside a function: calls from the global scope allocate)
function row_allocations(IK, i)
    n = length(IK.pdf)
    S = zeros(n, n)
    CI = CartesianIndices(IK.pdf.p)
    PIM.compute_row!(S, IK, i, CI[i], true)
    @allocated PIM.compute_row!(S, IK, i, CI[i], true)
end
@testset "The rows are computed without allocations" begin
    for axes in ((QuinticAxis(-4., 4., 21), QuinticAxis(-4., 4., 21)), (QuinticAxis(-4., 4., 21), ChebyshevAxis(-4., 4., 15)))
        IK = PathIntegration(duffing(), RK4(), 0.05, axes...; backtracing = ExplicitBacktracing(), extract_IK = Val(true))
        @test row_allocations(IK, 150) == 0
        K = PathIntegration(duffing(), RK4(), 0.05, axes...; backtracing = StrangSplitting(), extract_IK = Val(true))
        @test row_allocations(K.transport, 150) == 0
        @test all(IK -> row_allocations(IK, 150) == 0, K.diffusion)
    end
end

@testset "Back-tracing options" begin
    @test PathIntegration(SDE(f, g), RK4(), 0.05, CubicAxis(-3., 3., 21)).step_dynamics.steptracer isa PIM.SymbolicNewtonStepTracer # the default
    for bt in (NewtonBacktracing(), ExplicitBacktracing(), StrangSplitting())
        sdestep = SDEStep(SDE(f, g), RK4(), 0.05; backtracing = bt)
        @test PIM.backtracing(sdestep) == bt
        @test PIM.backtracing(PathIntegration(sdestep, 0.05, CubicAxis(-3., 3., 21); backtracing = bt).step_dynamics) == bt
    end
    # the back-tracing of an SDE step cannot be changed in PathIntegration
    @test_throws ArgumentError PathIntegration(SDEStep(SDE(f, g), RK4(), 0.05), 0.05, CubicAxis(-3., 3., 21); backtracing = ExplicitBacktracing())
    # not available for the vibro-impact oscillator
    sde_vio = SDE_VIO((x,p,t) -> -x[1], (x,p,t) -> 0.5, Wall(0.7, 0.0))
    @test_throws ArgumentError SDEStep(sde_vio, Euler(), 0.1; backtracing = ExplicitBacktracing())
    @test_throws ArgumentError SDEStep(sde_vio, Euler(), 0.1; backtracing = StrangSplitting())
end

end
