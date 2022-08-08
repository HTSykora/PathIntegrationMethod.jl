using PathIntegrationMethod
using Test

@testset "PathIntegrationMethod.jl" begin
    @time @testset "Interpolation" begin @test include("interpolationtest.jl") end
    @time @testset "Integration" begin @test include("integration_test.jl") end
    @time @testset "SDE type" begin @test include("sde_test.jl") end
    @time @testset "Steptracing test" begin @test include("steptracing_test.jl") end
    @time @testset "Scalar Cubic SDE" begin @test include("scalar_cubic_sde_test.jl") end


    # Write your tests here.
end