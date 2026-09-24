using PathIntegrationMethod
using Test

@testset "PathIntegrationMethod.jl" begin
    @time @testset "Interpolation" begin @test include("interpolationtest.jl") end
    @time @testset "Integration" begin @test include("integration_test.jl") end
    @time @testset "SDE type" begin @test include("sde_test.jl") end
    @time @testset "Steptracing test" begin @test include("steptracing_test.jl") end
    @time @testset "Scalar Cubic SDE" begin @test include("scalar_cubic_sde_test.jl") end
    @time @testset "Quadrature units" begin include("quadrature_unit_test.jl") end
    @time @testset "Interpolation units" begin include("interpolation_unit_test.jl") end
    @time @testset "Invariants" begin include("invariants_test.jl") end
    @time @testset "Time-periodic SDE" begin include("periodic_sde_test.jl") end
    @time @testset "Duffing oscillator (d=2)" begin include("duffing_2d_test.jl") end


    # Write your tests here.
end