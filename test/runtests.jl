using PathIntegrationMethod
using Test

@testset "PathIntegrationMethod.jl" begin
    @time @testset "Interpolation" begin @test include("interpolationtest.jl") end
    @time @testset "Integration" begin @test include("integration_test.jl") end
    @time @testset "SDE type" begin @test include("sde_test.jl") end
    @time @testset "Steptracing test" begin @test include("steptracing_test.jl") end
    @time @testset "Drift step" begin include("drift_step_test.jl") end
    @time @testset "Scalar Cubic SDE" begin @test include("scalar_cubic_sde_test.jl") end
    @time @testset "Quadrature units" begin include("quadrature_unit_test.jl") end
    @time @testset "Interpolation units" begin include("interpolation_unit_test.jl") end
    @time @testset "Invariants" begin include("invariants_test.jl") end
    @time @testset "Time-periodic SDE" begin include("periodic_sde_test.jl") end
    @time @testset "Duffing oscillator (d=2)" begin include("duffing_2d_test.jl") end
    @time @testset "Row kernels and sparse assembly" begin include("row_kernel_test.jl") end
    @time @testset "Steady state" begin include("steady_state_test.jl") end
    @time @testset "Parallel row computation" begin include("row_computation_test.jl") end
    @time @testset "Vibro-impact oscillator" begin include("vio_test.jl") end
    @time @testset "Type stability" begin include("type_stability_test.jl") end
    @time @testset "Back-tracing methods" begin include("backtracing_test.jl") end


    # Write your tests here.
end