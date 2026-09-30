module TypeStabilityTest
using PathIntegrationMethod, Test
const PIM = PathIntegrationMethod

# Duffing oscillator: dx = v dt, dv = (-2ζv + x - λx³) dt + σ dW
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]

# uniform and mixed axes: the axes (and the buffers of their interpolations) have different types in the mixed cases
axes_cases = (("quintic × quintic", () -> (QuinticAxis(-4., 4., 15), QuinticAxis(-4., 4., 15))),
              ("Chebyshev × Chebyshev", () -> (ChebyshevAxis(-4., 4., 11), ChebyshevAxis(-4., 4., 11))),
              ("quintic × Chebyshev", () -> (QuinticAxis(-4., 4., 15), ChebyshevAxis(-4., 4., 11))),
              ("Chebyshev × quintic", () -> (ChebyshevAxis(-4., 4., 11), QuinticAxis(-4., 4., 15))),
              ("cubic × quintic", () -> (CubicAxis(-4., 4., 15), QuinticAxis(-4., 4., 15))),
              ("Chebyshev × trigonometric", () -> (ChebyshevAxis(-4., 4., 11), TrigonometricAxis(-4., 4., 13))))

# (measured inside functions: calls from the global scope allocate)
function row_allocations(PI, i)
    CI = CartesianIndices(PI.pdf.p)
    S = PIM.storage_matrix(PI.stepMX[1]) # dense: the rows are written without buffers
    PIM.compute_row!(S, PI.IK, i, CI[i], true)
    @allocated PIM.compute_row!(S, PI.IK, i, CI[i], true)
end
evaluation_allocations(F, x::Vararg{Any,N}) where N = (F(x...); @allocated F(x...))

@testset "Type stability: $lbl" for (lbl, axes) in axes_cases
    PI = PathIntegration(SDE((f1, f2), g2, [0.5, 0.25, 1.0]), RK4(), 0.05, axes()...; stepMXtype = DenseMX(), rowcomputation = SerialRowComputation(), mPDF_IDs = [1, 2])
    # a row of the step matrix is computed without allocations
    @test row_allocations(PI, length(PI.pdf) ÷ 2) == 0
    advance!(PI)
    update_mPDFs!(PI)
    @test evaluation_allocations(PI, 0.1, 0.2) == 0
    @test evaluation_allocations(PI.marginal_pdfs[1], 0.1) == 0
    @test @inferred(integrate(PI.pdf)) ≈ 1
    @inferred integrate((x, v) -> x^2, PI.pdf)
    @test @inferred(integrate_diff(PI.pdf, PI.p_temp)) ≥ 0
end

end
