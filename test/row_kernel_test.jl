module RowKernelTest
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
# Linear system with d = 3
h1(x,p,t) = x[2]
h2(x,p,t) = -x[1] - x[2] + x[3]
h3(x,p,t) = -x[3] - x[2]
g3(x,p,t) = 0.8

# (label, SDE constructor, time stepping, time step, axes constructor)
sparse_systems = (("d=1 cubic", () -> SDE(f, g), RK2, 0.05, () -> (CubicAxis(-3., 3., 41),)),
                  ("d=2 quintic", () -> SDE((f1, f2), g2, copy(par)), Euler, 0.05, () -> (QuinticAxis(-4., 4., 15), QuinticAxis(-4., 4., 15))),
                  ("d=2 cubic × Chebyshev", () -> SDE((f1, f2), g2, copy(par)), Euler, 0.05, () -> (CubicAxis(-4., 4., 15), ChebyshevAxis(-4., 4., 11))),
                  ("d=3 cubic", () -> SDE((h1, h2, h3), g3), Euler, 0.1, () -> Tuple(CubicAxis(-3., 3., 9) for _ in 1:3)))
dense_systems = (("d=1 Chebyshev", () -> SDE(f, g), RK2, 0.05, () -> (ChebyshevAxis(-3., 3., 31),)),
                 ("d=1 trigonometric", () -> SDE(f, g), RK4, 0.05, () -> (TrigonometricAxis(-3., 3., 31),)),
                 ("d=2 Chebyshev", () -> SDE((f1, f2), g2, copy(par)), RK4, 0.05, () -> (ChebyshevAxis(-4., 4., 11), ChebyshevAxis(-4., 4., 11))),
                 ("d=3 Chebyshev", () -> SDE((h1, h2, h3), g3), Euler, 0.1, () -> Tuple(ChebyshevAxis(-3., 3., 7) for _ in 1:3)))

new_PI(sde, method, Δt, axes; kwargs...) = PathIntegration(sde(), method(), Δt, axes()...; kwargs...)
stepmatrix(args...; kwargs...) = Matrix(new_PI(args...; kwargs...).stepMX[1])

@testset "Sparse accumulator equals the generic row kernel: $lbl" for (lbl, sys...) in sparse_systems
    @test new_PI(sys...).IK.temp.kernel isa PIM.SparseAccumulator
    for stepMXtype in (SparseMX(), DenseMX())
        @test stepmatrix(sys...; stepMXtype) == stepmatrix(sys...; stepMXtype, generic_row_kernel = true)
    end
end

@testset "Dense tensor row kernel equals the generic row kernel: $lbl" for (lbl, sys...) in dense_systems
    @test new_PI(sys...).IK.temp.kernel isa PIM.DenseTensorKernel
    for stepMXtype in (SparseMX(), DenseMX())
        S = stepmatrix(sys...; stepMXtype)
        @test S ≈ stepmatrix(sys...; stepMXtype, generic_row_kernel = true) rtol = 1e-12
    end
end

@testset "The generic row kernel is used with QuadGK" begin
    PI = PathIntegration(SDE(f, g), RK4(), 0.05, CubicAxis(-3., 3., 21); discreteintegrator = QuadGKIntegrator())
    @test PI.IK.temp.kernel isa PIM.GenericRowKernel
end

@testset "Sparse step matrix equals the dense one without the small elements: $lbl" for (lbl, sys...) in sparse_systems[1:2]
    tol = 1e-6
    S_dense = stepmatrix(sys...; stepMXtype = DenseMX())
    @test stepmatrix(sys...; stepMXtype = SparseMX(sparse_tol = tol)) == map(x -> abs(x) > tol ? x : zero(x), S_dense)
end

@testset "Relative sparsity threshold" begin
    sys = sparse_systems[2][2:end]
    S0 = stepmatrix(sys...; stepMXtype = SparseMX(sparse_tol = 0.0))
    τ = 1e-8 * maximum(abs, S0)
    PI = new_PI(sys...; sparse_tol = 0.0, sparse_rtol = 1e-8)
    S_rel = Matrix(PI.stepMX[1])
    @test S_rel == map(x -> abs(x) > τ ? x : zero(x), S0)
    @test count(!iszero, S_rel) < count(!iszero, S0)
    recompute_stepMX!(PI)
    @test Matrix(PI.stepMX[1]) == S_rel
end

@testset "Int32 indices" begin
    sys = sparse_systems[2][2:end]
    PI64 = new_PI(sys...)
    PI32 = new_PI(sys...; index_type = Int32)
    @test eltype(PIM.rowvals(parent(PI32.stepMX[1]))) == Int32
    @test Matrix(PI32.stepMX[1]) == Matrix(PI64.stepMX[1])
    for _ in 1:3
        advance!(PI32)
        advance!(PI64)
    end
    @test PI32.pdf.p == PI64.pdf.p
    recompute_stepMX!(PI32)
    @test Matrix(PI32.stepMX[1]) == Matrix(PI64.stepMX[1])
    # too many nonzero elements for the index type
    @test_throws ArgumentError new_PI(sparse_systems[1][2:end]...; stepMXtype = SparseMX(index_type = Int8))
end

end
