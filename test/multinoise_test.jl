module MultiNoiseTest
using PathIntegrationMethod, Test, LinearAlgebra
const PIM = PathIntegrationMethod

# Linear systems with noise on several coordinates: the stationary PDF is the Gaussian N(0, P), A P + P Aᵀ + B Bᵀ = 0
quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
gaussian(P) = (x...) -> exp(-0.5*dot(collect(x), P \ collect(x)))
function stationary_error(PI, P)
    _, info = quiet(() -> steady_state!(PI))
    PI_exact = PathIntegration(PI.step_dynamics, PI.ts, PI.pdf.axes...; f_init = gaussian(P), pre_compute = false)
    PI_exact.pdf.p ./= integrate(PI_exact.pdf)
    integrate_diff(PI.pdf, PI_exact.pdf), info.residual
end

# d = 3, k = 2: dx = v dt, dv = (-x - 2ζv + y) dt + σ₁ dW₁, dy = -αy dt + σ₂ dW₂ (an oscillator driven by white noise and an Ornstein–Uhlenbeck process)
ζ, α, σ₁, σ₂ = 0.5, 1.0, 0.5, 0.7
h1(x,p,t) = x[2]
h2(x,p,t) = -x[1] - 2p[1]*x[2] + x[3]
h3(x,p,t) = -p[2]*x[3]
s2(x,p,t) = p[3]
s3(x,p,t) = p[4]
sde3(par = [ζ, α, σ₁, σ₂]) = SDE((h1, h2, h3), (s2, s3), par)
P3 = lyap([0 1 0; -1 -2ζ 1; 0 0 -α], Diagonal([0, σ₁^2, σ₂^2]))
axes3(N) = Tuple(QuinticAxis(-L, L, N) for L in 4.5 .* sqrt.(diag(P3)))
# d = 2, k = 1: dx = -x dt + σ₁ dW₁, dy = (x - y) dt + σ₂ dW₂
u1(x,p,t) = -x[1]
u2(x,p,t) = x[1] - x[2]
sde2 = SDE((u1, u2), ((x,p,t) -> σ₁, (x,p,t) -> σ₂))
P2 = lyap([-1.0 0; 1 -1], Diagonal([σ₁^2, σ₂^2]))
axes2(N) = Tuple(QuinticAxis(-L, L, N) for L in 4.5 .* sqrt.(diag(P2)))
stepmatrices(PI) = [Matrix(S) for S in PI.stepMX]

@testset "Diagonal noise on several coordinates" begin
    @test PIM.get_dkm(sde3()) == (3, 2, 2)
    @test PIM.get_dkm(sde2) == (2, 1, 2)
    PI = PathIntegration(sde3(), RK4(), 0.05, axes3(17)...)
    @test PI.IK.discreteintegrator isa PIM.GaussianDiscreteIntegrator{2}
    @test length(PI.IK.discreteintegrator.x) == 7^2 # the default: Gauss–Hermite with 7 nodes in each noisy coordinate
end

@testset "Exact stationary PDF: d = $d, k = $k, $bt" for (d, k, new_PI, P) in ((3, 2, bt -> PathIntegration(sde3(), RK4(), 0.05, axes3(17)...; backtracing = bt), P3),
                                                                                (2, 1, bt -> PathIntegration(sde2, RK4(), 0.05, axes2(31)...; backtracing = bt), P2)),
                                                                          bt in (NewtonBacktracing(), ExplicitBacktracing())
    e, residual = stationary_error(new_PI(bt), P)
    @test residual < 1e-10
    @test e < 0.05 # the error of the time step (O(Δt))
end

@testset "First order convergence in the time step" begin
    e1, _ = stationary_error(PathIntegration(sde3(), RK4(), 0.1, axes3(17)...; backtracing = ExplicitBacktracing()), P3)
    e2, _ = stationary_error(PathIntegration(sde3(), RK4(), 0.05, axes3(17)...; backtracing = ExplicitBacktracing()), P3)
    @test 1.6 < e1 / e2 < 2.4
end

@testset "Quadrature rules" begin
    new_PI(di) = PathIntegration(sde3(), RK4(), 0.05, axes3(17)...; backtracing = ExplicitBacktracing(), discreteintegrator = di)
    PI_GL = new_PI(GaussLegendreIntegrator(21, dim = 2))
    quiet(() -> steady_state!(PI_GL))
    for di in (GaussHermiteIntegrator(5, dim = 2), GaussHermiteIntegrator(7, dim = 2))
        PI = new_PI(di)
        quiet(() -> steady_state!(PI))
        @test integrate_diff(PI.pdf, PI_GL.pdf.p) < 6e-3 # (the error of the time step is about 4e-2)
    end
    # a one-dimensional rule is used in every noisy coordinate
    @test stepmatrices(new_PI(GaussLegendreIntegrator(9))) == stepmatrices(new_PI(GaussLegendreIntegrator(9, dim = 2)))
    @test length(new_PI(GaussHermiteIntegrator((5, 7))).IK.discreteintegrator.x) == 35
end

@testset "Quadrature nodes" begin
    # Gauss–Hermite: exact for a Gaussian times a polynomial of degree 2N - 1
    gh = PIM.DiscreteIntegrator(GaussHermiteIntegrator(4), zeros(0), CubicAxis(-1., 1., 5))
    c, s = 0.3, 0.2
    PIM.rescale_to_gaussian!(gh, (c,), (s,))
    f(z) = exp(-(z - c)^2 / (2s^2)) / sqrt(2π*s^2) * (1 + z + z^2 + z^7)
    exact = 1 + c + (c^2 + s^2) + (c^7 + 21c^5*s^2 + 105c^3*s^4 + 105c*s^6)
    @test sum(w * f(z) for (w, z) in zip(gh.w, gh.x)) ≈ exact
    # the tensor product of the one-dimensional rules (rescaled in the same way: the outer nodes are mapped to the ends of the window)
    gl = PIM.DiscreteIntegrator(GaussLegendreIntegrator(3, dim = 2), zeros(0), CubicAxis(-1., 1., 5), CubicAxis(-1., 1., 5))
    PIM.rescale_to_limits!(gl, (0.0, -1.0), (2.0, 3.0))
    @test gl.x[1] isa PIM.SVector{2}
    gl1(start, stop) = (q = PIM.DiscreteIntegrator(GaussLegendreIntegrator(3), zeros(0), CubicAxis(-1., 1., 5)); PIM.rescale_to_limits!(q, start, stop); q)
    q1, q2 = gl1(0.0, 2.0), gl1(-1.0, 3.0)
    @test sum(w * x[1]^5 * x[2]^3 for (w, x) in zip(gl.w, gl.x)) ≈ sum(q1.w .* q1.x .^ 5) * sum(q2.w .* q2.x .^ 3)
    @test (q1.x[1], q1.x[end]) == (0.0, 2.0)
end

@testset "Row computations: $bt" for bt in (NewtonBacktracing(), ExplicitBacktracing())
    new_PI(rc) = PathIntegration(sde3(), RK4(), 0.05, CubicAxis(-2., 2., 11), CubicAxis(-2., 2., 11), CubicAxis(-2., 2., 11); backtracing = bt, rowcomputation = rc)
    S = stepmatrices(new_PI(SerialRowComputation()))
    @test stepmatrices(new_PI(ThreadedRowComputation(3))) == S
    @test stepmatrices(new_PI(BatchRowComputation(3))) == S
    # recompute with new parameters
    PI = new_PI(ThreadedRowComputation(3))
    recompute_PI!(PI; par = [0.4, 0.8, 0.5, 0.6])
    @test stepmatrices(PI) == stepmatrices(PathIntegration(sde3([0.4, 0.8, 0.5, 0.6]), RK4(), 0.05, CubicAxis(-2., 2., 11), CubicAxis(-2., 2., 11), CubicAxis(-2., 2., 11);
        backtracing = bt, rowcomputation = SerialRowComputation()))
end

# (measured inside a function: calls from the global scope allocate)
function row_allocations(IK, i)
    n = length(IK.pdf)
    S = zeros(n, n)
    CI = CartesianIndices(IK.pdf.p)
    PIM.compute_row!(S, IK, i, CI[i], true)
    @allocated PIM.compute_row!(S, IK, i, CI[i], true)
end
@testset "The rows are computed without allocations: $bt" for bt in (NewtonBacktracing(), ExplicitBacktracing())
    IK = PathIntegration(sde3(), RK4(), 0.05, CubicAxis(-2., 2., 11), CubicAxis(-2., 2., 11), CubicAxis(-2., 2., 11); backtracing = bt, extract_IK = Val(true))
    @test row_allocations(IK, 600) == 0
end

@testset "Gauss–Hermite with noise on the last coordinate" begin
    f1(x,p,t) = x[2]
    f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
    g2(x,p,t) = p[3]
    new_PI(di) = PathIntegration(SDE((f1, f2), g2, [0.5, 0.25, 1.0]), RK4(), 0.05, QuinticAxis(-4., 4., 41), QuinticAxis(-4., 4., 41); discreteintegrator = di)
    PI_GL, PI_GH = new_PI(GaussLegendreIntegrator(31)), new_PI(GaussHermiteIntegrator(7))
    quiet(() -> steady_state!(PI_GL)); quiet(() -> steady_state!(PI_GH))
    @test integrate_diff(PI_GH.pdf, PI_GL.pdf.p) < 1e-3
end

@testset "StrangSplitting" begin
    # the diffusion of each noisy coordinate, one-dimensional integrals (Gauss–Legendre by default)
    K = PathIntegration(sde3(), RK4(), 0.05, axes3(9)...; backtracing = StrangSplitting(), extract_IK = Val(true))
    @test length(K.diffusion) == 2
    @test all(IK -> IK.discreteintegrator isa PIM.IntervalDiscreteIntegrator{1} && length(IK.discreteintegrator.x) == 31, K.diffusion)
    # a rule for two noisy coordinates gives the rule of each coordinate
    K = PathIntegration(sde3(), RK4(), 0.05, axes3(9)...; backtracing = StrangSplitting(), discreteintegrator = GaussLegendreIntegrator((15, 21)), extract_IK = Val(true))
    @test map(IK -> length(IK.discreteintegrator.x), K.diffusion) == (15, 21)

    # second order for additive noise (d = 2, noise on both coordinates)
    new_PI(Δt) = PathIntegration(sde2, RK4(), Δt, axes2(41)...; backtracing = StrangSplitting())
    e1, residual = stationary_error(new_PI(0.1), P2)
    e2, _ = stationary_error(new_PI(0.05), P2)
    @test residual < 1e-10
    @test e1 / e2 > 3
    @test e2 < 1e-3
    # d = 3: much more accurate than the Maruyama approximation of the whole time step (limited by the interpolation on this coarse grid)
    e_strang, _ = stationary_error(PathIntegration(sde3(), RK4(), 0.1, axes3(21)...; backtracing = StrangSplitting()), P3)
    e_explicit, _ = stationary_error(PathIntegration(sde3(), RK4(), 0.1, axes3(21)...; backtracing = ExplicitBacktracing()), P3)
    @test e_strang < e_explicit / 5

    # row computations and recomputation
    new_PI3(par, rc) = PathIntegration(sde3(par), RK4(), 0.05, CubicAxis(-2., 2., 9), CubicAxis(-2., 2., 9), CubicAxis(-2., 2., 9); backtracing = StrangSplitting(), rowcomputation = rc)
    par = [ζ, α, σ₁, σ₂]
    S = stepmatrices(new_PI3(copy(par), SerialRowComputation()))
    @test stepmatrices(new_PI3(copy(par), ThreadedRowComputation(3))) == S
    @test stepmatrices(new_PI3(copy(par), BatchRowComputation(3))) == S
    PI = new_PI3(copy(par), ThreadedRowComputation(3))
    recompute_PI!(PI; par = 0.9 .* par)
    @test stepmatrices(PI) == stepmatrices(new_PI3(0.9 .* par, SerialRowComputation()))
    # the rows of the transport and of the diffusion steps are computed without allocations
    K = PathIntegration(sde3(), RK4(), 0.05, CubicAxis(-2., 2., 11), CubicAxis(-2., 2., 11), CubicAxis(-2., 2., 11); backtracing = StrangSplitting(), extract_IK = Val(true))
    @test row_allocations(K.transport, 600) == 0
    @test all(IK -> row_allocations(IK, 600) == 0, K.diffusion)
end

@testset "Not available with noise on several coordinates" begin
    @test_throws ArgumentError PathIntegration(sde3(), RK4(), 0.05, axes3(9)...; discreteintegrator = QuadGKIntegrator())
    @test_throws ArgumentError PathIntegration(sde3(), RK4(), 0.05, axes3(9)...; discreteintegrator = QuadGKIntegrator(), backtracing = StrangSplitting())
    @test_throws ArgumentError PathIntegration(sde3(), RK4(), 0.05, axes3(9)...; discreteintegrator = GaussLegendreIntegrator((5, 5, 5)))
    @test_throws ArgumentError PathIntegration(sde3(), RK4(), 0.05, axes3(9)...; discreteintegrator = GaussLegendreIntegrator((5, 5, 5)), backtracing = StrangSplitting())
    @test_throws ArgumentError PathIntegration(sde3(), RK4(), 0.05, axes3(9)...; smart_integration = false)
end

end
