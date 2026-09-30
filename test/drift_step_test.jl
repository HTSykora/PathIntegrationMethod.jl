module DriftStepTest
using PathIntegrationMethod, Test
const PIM = PathIntegrationMethod

# Time-dependent drift: dx = (p₁ + p₂t + p₃t² + p₄t³) dt + σ dW
ft(x,p,t) = p[1] + p[2]*t + p[3]*t^2 + p[4]*t^3
g(x,p,t) = 0.5
# The drift step is x₁ = x₀ + ∫f(t)dt over [t₀, t₁]: RK2 integrates it exactly up to linear, RK4 up to cubic f(t)
F(t, p) = p[1]*t + p[2]*t^2/2 + p[3]*t^3/3 + p[4]*t^4/4
t0, t1, x0 = 1.0, 1.1, 0.3

# (label, time stepping, parameters)
methods = (("Euler", Euler(), [0.5, 0.0, 0.0, 0.0]),
           ("RK2(α = 1//2)", RK2(α = 1//2), [0.5, 2.0, 0.0, 0.0]),
           ("RK2(α = 2//3)", RK2(α = 2//3), [0.5, 2.0, 0.0, 0.0]),
           ("RK2(α = 1)", RK2(α = 1), [0.5, 2.0, 0.0, 0.0]),
           ("RK4", RK4(), [0.5, 2.0, -1.0, 3.0]))

@testset "Stage times of the drift step: $lbl" for (lbl, method, par) in methods
    sde = SDE(ft, g, par)
    x_exact = x0 + F(t1, par) - F(t0, par)
    step = SDEStep(sde, method, [x0], [0.0], t0, t1)
    PIM.eval_driftstep!(step)
    @test step.x1[1] ≈ x_exact
    # the symbolic drift step (compiled for the Newton back-tracing)
    @test PIM.eval_driftstep_xI_sym(sde, step.method, [x0], par, t0, t1)[1] ≈ x_exact
end

@testset "Mean of the PDF after a step: $lbl" for (lbl, method, par) in methods
    # the Gaussian initial PDF is shifted by the drift step
    x_exact = x0 + F(t1, par) - F(t0, par)
    PI = PathIntegration(SDE(ft, g, par), method, [t0, t1], ChebyshevAxis(-2., 4., 61); f_init = x -> exp(-(x - x0)^2 / (2*0.2^2)))
    advance!(PI)
    ax = PI.pdf.axes[1]
    @test sum(ax.wts .* ax.xs .* PI.pdf.p) / sum(ax.wts .* PI.pdf.p) ≈ x_exact atol = 1e-8
end

end
