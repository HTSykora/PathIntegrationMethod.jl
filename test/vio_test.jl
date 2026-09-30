module VIOTest
using PathIntegrationMethod, Test

# Linear oscillator with a rigid wall at its equilibrium x = 0: dx = v dt, dv = (-2ζv - x) dt + σ dW, x ≥ 0
f2(x,p,t) = -2p[1]*x[2] - x[1]
g2(x,p,t) = p[2]
ζ, σ = 0.25, 0.5
new_PI(wall; kwargs...) = PathIntegration(SDE_VIO(f2, g2, wall, [ζ, σ]), RK4(), 0.05, QuinticAxis(0., 3., 41), QuinticAxis(-3., 3., 61); kwargs...)
stepmatrices(PI) = [Matrix(S) for S in PI.stepMX]

@testset "Elastic impacts: exact stationary PDF" begin
    # with r = 1 the stationary PDF is the Gaussian of the oscillator without the wall, restricted to x ≥ 0
    s² = σ^2 / (4ζ)
    p_exact(x, v) = 2 * exp(-(x^2 + v^2) / (2s²)) / (2π*s²)
    PI = new_PI(Wall(1.0, 0.0))
    PI_exact = new_PI(Wall(1.0, 0.0); f_init = p_exact, pre_compute = false)
    _, info = steady_state!(PI)
    @test info.residual < 1e-10
    @test integrate(PI.pdf) ≈ 1
    @test integrate_diff(PI.pdf, PI_exact.pdf) < 0.02 # the error is O(Δt): 0.012

    # the step matrix does not depend on the row computation
    S_serial = stepmatrices(PI)
    for rowcomputation in (ThreadedRowComputation(3), BatchRowComputation(3))
        recompute_stepMX!(PI; rowcomputation)
        @test stepmatrices(PI) == S_serial
    end
end

@testset "Restitution coefficient as a function of the impact velocity" begin
    # a constant function gives the same step matrix as the constant coefficient
    @test stepmatrices(new_PI(Wall(v -> 0.7 + 0*v, 0.0))) == stepmatrices(new_PI(Wall(0.7, 0.0)))

    PI = new_PI(Wall(v -> 0.9 - 0.2tanh(v), 0.0))
    _, info = steady_state!(PI)
    @test info.residual < 1e-10
    @test integrate(PI.pdf) ≈ 1
    @test minimum(PI.pdf.p) > -1e-8
end

end
