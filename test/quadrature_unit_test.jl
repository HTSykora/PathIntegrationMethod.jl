module QuadratureUnitTest
using PathIntegrationMethod, Test
const PIM = PathIntegrationMethod

a, b = -1.0, 2.0
exact(k) = (b^(k+1) - a^(k+1)) / (k+1) # ∫ₐᵇ xᵏ dx

@testset "Newton–Cotes weights" begin
    # Order 1 is exact up to degree 1, orders 2 and 3 up to degree 3
    # l = 11:16 covers every end correction (remainder mod(l-1, order))
    for order in 1:3, l in 11:16
        w = collect(PIM.NewtonCotesWeights(order, l, (b - a)/(l - 1)))
        x = LinRange(a, b, l)
        err = maximum(abs(sum(w .* x.^k) - exact(k)) for k in 0:(order == 1 ? 1 : 3))
        @test err < 1e-12
    end
end

@testset "Clenshaw–Curtis weights" begin
    # exact up to degree n-1
    for n in (8, 9)
        x = PIM.chebygrid(Float64, a, b, n)
        w = PIM.clenshawcurtisweights(Float64, a, b, n)
        @test maximum(abs(sum(w .* x.^k) - exact(k)) for k in 0:n-1) < 1e-12
    end
end

@testset "QuadGKIntegrator" begin
    # Sparse interpolations only write a few elements of the integrand
    f!(v, x) = (v[2] = exp(-x^2); v)
    I = quadgk(x -> exp(-x^2), a, b)[1]
    di = DiscreteIntegrator(QuadGKIntegrator(), zeros(5), CubicAxis(a, b, 11))
    di(f!)
    @test di.res[2] ≈ I
    @test all(iszero, di.res[[1, 3, 4, 5]])

    res = ones(5)
    di(f!, res; Q_reinit_res = false) # accumulate
    @test res[2] ≈ 1 + I
    @test all(==(1), res[[1, 3, 4, 5]])

    PIM.rescale_to_limits!(di, 1.0, 1.0) # empty interval
    di(f!)
    @test all(iszero, di.res)
end

@testset "Rescaling of the quadrature nodes does not depend on the previous rescalings" begin
    ax = CubicAxis(a, b, 11)
    di = DiscreteIntegrator(GaussLegendreIntegrator(31), zeros(11), ax)
    di_fresh = DiscreteIntegrator(GaussLegendreIntegrator(31), zeros(11), ax)
    for (start, stop) in ((0.3, 0.30001), (-0.7, 1.9), (1.2, 1.7))
        PIM.rescale_to_limits!(di, start, stop)
    end
    PIM.rescale_to_limits!(di_fresh, 1.2, 1.7)
    @test di.x == di_fresh.x
    @test di.w == di_fresh.w
    @test di.x[1] == 1.2 && di.x[end] ≈ 1.7
end

@testset "Discrete integrators in PathIntegration" begin
    # Newton–Cotes type integrators store equidistant nodes (they used to fail with smart integration)
    f(x,p,t) = x[1] - x[1]^3
    g(x,p,t) = sqrt(2)
    S_GL = Matrix(PathIntegration(SDE(f, g), RK4(), 0.01, CubicAxis(-3., 3., 41)).stepMX[1])
    for di in (NewtonCotesIntegrator(31, 2), TrapezoidalIntegrator(31), ClenshawCurtisIntegrator(31))
        S = Matrix(PathIntegration(SDE(f, g), RK4(), 0.01, CubicAxis(-3., 3., 41); discreteintegrator = di).stepMX[1])
        @test maximum(abs, S - S_GL) < 0.02 * maximum(abs, S_GL)
    end
end

@testset "QuadGKIntegrator in PathIntegration" begin
    f(x,p,t) = x[1] - x[1]^3
    g(x,p,t) = sqrt(2)
    for axis in (CubicAxis, ChebyshevAxis)
        S_GL = Matrix(PathIntegration(SDE(f, g), RK4(), 0.01, axis(-3., 3., 41)).stepMX[1])
        S_QGK = Matrix(PathIntegration(SDE(f, g), RK4(), 0.01, axis(-3., 3., 41); discreteintegrator = QuadGKIntegrator()).stepMX[1])
        @test all(isfinite, S_QGK)
        @test maximum(abs, S_QGK - S_GL) < 0.02 * maximum(abs, S_GL)
    end
end

end
