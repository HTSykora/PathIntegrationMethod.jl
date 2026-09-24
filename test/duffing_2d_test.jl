module Duffing2DTest
using PathIntegrationMethod, Test
const PIM = PathIntegrationMethod

# Duffing oscillator (Eq. 31 in Sykora, Kuske & Yurchenko, 2022)
#   dx = v dt,   dv = (-2ζv + x - λx³) dt + σ dW
# with the steady-state PDF p(x,v) ∝ exp(C₀(x²/2 - λx⁴/4 - v²/2)), C₀ = 4ζ/σ²
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]
par = [0.5, 0.25, 1.0] # ζ, λ, σ

C₀ = 4par[1]/par[3]^2
xmin, xmax = -4.0, 4.0
_p_x(x) = exp(C₀*(x^2/2 - par[2]*x^4/4))
_p_v(v) = exp(-C₀*v^2/2)
Z_x = quadgk(_p_x, xmin, xmax)[1]
Z_v = quadgk(_p_v, xmin, xmax)[1]
p_x(x) = _p_x(x) / Z_x # marginal PDFs normalised on the computational domain
p_v(v) = _p_v(v) / Z_v

xs = LinRange(xmin, xmax, 201)
err_2D(PI) = sum(abs(PI.pdf(x, v) - p_x(x)*p_v(v)) for x in xs, v in xs) * step(xs)^2
err_1D(f, p) = sum(abs(f(x) - p(x)) for x in xs) * step(xs)

err_margin = 0.2 # allowed relative increase of ε₁ compared to the reference errors
Δt, Tmax = 0.02, 40.0

@testset "Steady state vs analytic PDF" begin
    # (axis, N, time stepping, reference ε₁ of the joint PDF, marginal PDF ID => reference ε₁)
    cases = ((QuinticAxis, 41, Euler, 0.0346, [2 => 0.0202]),
             (QuinticAxis, 41, RK4, 0.0101, [1 => 0.00182, 2 => 0.00978]),
             (ChebyshevAxis, 31, RK4, 0.00969, []))
    for (axis, N, method, err_ref, mPDF_refs) in cases
        # a single marginal PDF is stored on its own, multiple ones in a tuple
        mPDF_IDs = isempty(mPDF_refs) ? nothing : length(mPDF_refs) == 1 ? first(mPDF_refs[1]) : first.(mPDF_refs)
        sde = SDE((f1, f2), g2, copy(par))
        PI = PathIntegration(sde, method(), Δt, axis(xmin, xmax, N), axis(xmin, xmax, N); mPDF_IDs = mPDF_IDs)
        advance_till_converged!(PI; Tmax = Tmax)
        @test integrate(PI.pdf) ≈ 1
        @test err_2D(PI) < (1 + err_margin) * err_ref

        isempty(mPDF_refs) && continue
        # Marginal PDFs
        update_mPDFs!(PI)
        mpdfs = PI.marginal_pdfs isa Tuple ? PI.marginal_pdfs : (PI.marginal_pdfs,)
        for (mpdf, (ID, mpdf_err_ref)) in zip(mpdfs, mPDF_refs)
            other = 3 - ID
            w = PI.pdf.axes[other].wts
            p_direct = [sum(w[j] * (ID == 1 ? PI.pdf.p[i, j] : PI.pdf.p[j, i]) for j in eachindex(w)) for i in eachindex(PI.pdf.axes[ID])]
            @test mpdf.pdf.p ≈ p_direct
            @test integrate(mpdf.pdf) ≈ 1
            @test err_1D(mpdf, ID == 1 ? p_x : p_v) < (1 + err_margin) * mpdf_err_ref
        end
    end
end

@testset "Jacobian correction det J" begin
    # 1/|det J_I| of the drift step x_I ↦ η_I(x_I, x_II) against central finite differences
    h = 1e-6
    x0 = [0.7, -0.4]
    for method in (Euler, RK2, RK4)
        st = SDEStep(SDE((f1, f2), g2, copy(par)), method(), zeros(2), zeros(2), 0.0, 0.1)
        η_I(x_I) = (st.x0 .= (x_I, x0[2]); PIM.eval_driftstep!(st); st.x1[1])
        J = (η_I(x0[1] + h) - η_I(x0[1] - h)) / 2h
        st.x0 .= x0
        @test PIM.get_detJinv(st) ≈ 1/abs(J) rtol = 1e-7
    end

    # d = 3 with two states without noise: J_I is 2×2
    h1(x,p,t) = x[2]^2 + x[1]^2 + x[3]^2 + t
    h2(x,p,t) = -x[1] - p[2]*x[1]^3 - 2p[1]*x[2] + sin(t) + p[2]*x[3]^3 + 2p[1]*x[3]
    h3(x,p,t) = -x[1] - p[2]*x[1]^3 - 2p[1]*x[2] + cos(t) - p[2]*x[3]^2
    g3(x,p,t) = sqrt(2)
    x0 = [0.3, -0.2, 0.5]
    for method in (Euler, RK2)
        st = SDEStep(SDE((h1, h2, h3), g3, [0.1, 0.1]), method(), zeros(3), zeros(3), 0.6, 0.7)
        η_I(x_I) = (st.x0 .= (x_I[1], x_I[2], x0[3]); PIM.eval_driftstep!(st); st.x1[1:2])
        J = hcat(((η_I(x0[1:2] .+ h .* e) .- η_I(x0[1:2] .- h .* e)) ./ 2h for e in ([1.0, 0.0], [0.0, 1.0]))...)
        st.x0 .= x0
        @test PIM.get_detJinv(st) ≈ 1/abs(det(J)) rtol = 1e-7
    end
end

end
