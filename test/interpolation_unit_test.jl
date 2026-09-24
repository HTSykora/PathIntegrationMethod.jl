module InterpolationUnitTest
using PathIntegrationMethod, Test
const PIM = PathIntegrationMethod

axistypes = (LinearAxis, CubicAxis, QuinticAxis, ChebyshevAxis, TrigonometricAxis)

@testset "Node values are reproduced: $axis" for axis in axistypes
    a1, a2 = axis(-1., 2., 12), axis(0., 1., 9)
    F1 = InterpolatedFunction(a1)
    F1.p .= rand(size(F1.p)...)
    @test all(F1(a1[i]) ≈ F1.p[i] for i in eachindex(a1))

    F2 = InterpolatedFunction(a1, a2)
    F2.p .= rand(size(F2.p)...)
    @test all(F2(a1[i], a2[j]) ≈ F2.p[i, j] for i in eachindex(a1), j in eachindex(a2))

    a3 = (axis === ChebyshevAxis ? CubicAxis : ChebyshevAxis)(0., 1., 9) # mixed interpolation types
    F3 = InterpolatedFunction(a1, a3)
    F3.p .= rand(size(F3.p)...)
    @test all(F3(a1[i], a3[j]) ≈ F3.p[i, j] for i in eachindex(a1), j in eachindex(a3))
end

@testset "Trigonometric interpolation, N = $N" for N in (16, 17)
    axis = TrigonometricAxis(-1., 2., N)
    L = N * (axis[2] - axis[1]) # period of the interpolation
    h(x) = 1 + sin(2π*(x + 1)/L) + 0.5cos(4π*(x + 1)/L)
    F = InterpolatedFunction(axis; f = h)
    @test maximum(abs(F(x) - h(x)) for x in LinRange(-1., 2., 301)) < 1e-12
end

# Basis function values and the interpolated value from them
function basis(axis, x; kwargs...)
    vals = axis.temp isa AbstractVector ? similar(axis.temp) : zero(axis.temp)
    PIM.basefun_vals_safe!(vals, axis, x; kwargs...)
    vals
end
itp_value(vals::AbstractVector, p) = sum(vals .* p)
itp_value(vals, p) = sum(p[i] * v for (i, v) in zip(vals.idxs, vals.val) if !iszero(v); init = zero(eltype(p)))

@testset "Extrapolation flags: $axis" for axis in axistypes
    lin(x) = 2x + 1
    a = axis(-1., 2., 12)
    p = lin.(a)
    Δ = (a[end] - a[1]) / 11
    below, above = a[1] - 0.3Δ, a[end] + 0.3Δ

    # default: zero outside the grid
    @test itp_value(basis(a, below), p) == 0
    @test itp_value(basis(a, above), p) == 0
    # zero_extrapolation = false: the boundary values are continued
    @test itp_value(basis(a, below; zero_extrapolation = false), p) ≈ p[1]
    @test itp_value(basis(a, above; zero_extrapolation = false), p) ≈ p[end]
    # allow_extrapolation = true: the boundary interpolation is continued (exact for a linear function)
    axis === TrigonometricAxis && continue # periodic interpolation
    @test itp_value(basis(a, below; allow_extrapolation = true), p) ≈ lin(below)
    if axis === ChebyshevAxis
        @test itp_value(basis(a, above; allow_extrapolation = true), p) ≈ lin(above)
    else
        # sparse interpolations: the stencil indices point beyond the last node
        @test_broken itp_value(basis(a, above; allow_extrapolation = true), p) ≈ lin(above)
    end
end

@testset "Extrapolation flags in InterpolatedFunction evaluation" begin
    F = InterpolatedFunction(CubicAxis(-1., 2., 12); f = x -> 2x + 1)
    @test F(2.1) == 0
    # the keywords are not forwarded by `interpolate`: always zero extrapolation
    @test_broken F(2.1; zero_extrapolation = false) ≈ F.p[end]
end

end
