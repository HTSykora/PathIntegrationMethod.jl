module SteadyStateTest
using PathIntegrationMethod, Test
const PIM = PathIntegrationMethod

# Scalar system: dx = (x - x³) dt + √2 dW
f(x,p,t) = x[1] - x[1]^3
g(x,p,t) = sqrt(2)
# Duffing oscillator: dx = v dt, dv = (-2ζv + x - λx³) dt + σ dW
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]
par = [0.5, 0.25, 1.0]
# Scalar system with periodic forcing: dx = (x - x³ + A sin(2πt/T)) dt + √2 dW
fp(x,p,t) = x[1] - x[1]^3 + p[1]*sin(2π*t/p[2])

# (label, PathIntegration constructor, methods to compare)
systems = (("d=1 cubic", () -> PathIntegration(SDE(f, g), RK4(), 0.01, CubicAxis(-3., 3., 71)), (:auto, :eigen, :lu, :arnoldi)),
           ("d=1 Chebyshev", () -> PathIntegration(SDE(f, g), RK4(), 0.01, ChebyshevAxis(-3., 3., 71)), (:auto, :eigen, :lu, :arnoldi)),
           ("d=2 quintic", () -> PathIntegration(SDE((f1, f2), g2, copy(par)), RK4(), 0.02, QuinticAxis(-4., 4., 41), QuinticAxis(-4., 4., 41)), (:auto, :lu, :arnoldi)),
           ("periodic", () -> PathIntegration(SDE(fp, g, [1.0, 0.5]), Euler(), collect(range(0, 0.5, length = 6)), CubicAxis(-3., 3., 41)), (:auto, :eigen, :arnoldi)),
           ("periodic Chebyshev", () -> PathIntegration(SDE(fp, g, [1.0, 0.5]), RK4(), collect(range(0, 0.5, length = 11)), ChebyshevAxis(-3., 3., 61)), (:auto, :eigen, :arnoldi)))

@testset "Steady state: $lbl" for (lbl, new_PI, methods) in systems
    PI, info = @test_logs steady_state!(new_PI()) # no warnings
    p = copy(PI.pdf.p)
    @test integrate(PI.pdf) ≈ 1
    @test info.residual < 1e-9
    @test abs(info.λ - 1) < 0.5
    @test PI.step_idx == 0
    @test PI.t == 0

    # p is a fixed point of advance! (for a whole period in the time-periodic case)
    for _ in eachindex(PI.stepMX)
        advance!(PI)
    end
    @test integrate_diff(PI.pdf, p) < 1e-9

    # the methods give the same PDF
    for method in methods
        PI_m, info_m = steady_state!(new_PI(); method)
        @test integrate_diff(PI_m.pdf, p) < 1e-8
        @test info_m.λ ≈ info.λ
        @test info_m.λ₂ ≈ info.λ₂
    end
end

@testset "Second eigenvalue" begin
    PI = systems[1][2]()
    λs = eigvals(Matrix(PI.stepMX[1]))
    λ₂ = sort(λs, by = μ -> abs(μ - 1))[2]
    for method in systems[1][3]
        @test steady_state!(systems[1][2](); method)[2].λ₂ ≈ λ₂
    end
end

@testset "Warning for several eigenvalues near 1" begin
    # bistable system with weak noise: the rare transitions between the wells give a second eigenvalue near 1
    new_PI() = PathIntegration(SDE(f, (x,p,t) -> 0.35), RK4(), 0.01, CubicAxis(-2., 2., 201))
    for method in (:eigen, :lu, :arnoldi)
        _, info = @test_logs (:warn, r"Several eigenvalues are near 1") steady_state!(new_PI(); method)
        @test abs(info.λ₂ - 1) < 1e-4
    end
    @test_logs steady_state!(new_PI(); near_one_tol = 1e-6) # a stricter definition of "near"
end

@testset "Steady state equals the converged advance!" begin
    PI_ss, _ = steady_state!(systems[1][2]())
    PI_adv = systems[1][2]()
    advance_till_converged!(PI_adv; rtol = 1e-10, Tmax = 100.0)
    @test integrate_diff(PI_ss.pdf, PI_adv.pdf) < 1e-8
end

@testset "Unsupported methods" begin
    @test_throws ArgumentError steady_state!(systems[4][2](); method = :lu) # several step matrices
    @test_throws ArgumentError steady_state!(systems[1][2](); method = :power)
end

end
