module InvariantsTest
using PathIntegrationMethod, Test
const PIM = PathIntegrationMethod

# Scalar system: dx = (x - p₁x³) dt + √2 dW
f(x,p,t) = x[1] - p[1]*x[1]^3
g(x,p,t) = sqrt(2)
# Duffing oscillator: dx = v dt, dv = (-2ζv + x - λx³) dt + σ dW
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]

# (label, SDE constructor with parameters par, parameters, axes constructor, time stepping, tolerance of the probability lost in a step)
systems = (("d=1 cubic", par -> SDE(f, g, par), [1.0], () -> (CubicAxis(-3., 3., 41),), RK2, 1e-4),
           ("d=1 Chebyshev", par -> SDE(f, g, par), [1.0], () -> (ChebyshevAxis(-3., 3., 31),), RK2, 1e-4),
           ("d=2 quintic", par -> SDE((f1, f2), g2, par), [0.5, 0.25, 1.0], () -> (QuinticAxis(-4., 4., 15), QuinticAxis(-4., 4., 15)), Euler, 2e-3))
Δt = 0.05

@testset "Step matrix representation: $lbl" for (lbl, sde, par, axes, method, mass_tol) in systems
    S_threaded = Matrix(PathIntegration(sde(copy(par)), method(), Δt, axes()...; stepMXtype = SparseMX(threaded = true)).stepMX[1])
    S_serial = Matrix(PathIntegration(sde(copy(par)), method(), Δt, axes()...; stepMXtype = SparseMX(threaded = false)).stepMX[1])
    S_dense = Matrix(PathIntegration(sde(copy(par)), method(), Δt, axes()...; stepMXtype = DenseMX()).stepMX[1])
    @test S_threaded == S_serial
    @test maximum(abs, S_dense - S_threaded) ≤ 1e-6 # default sparse_tol
    @test maximum(abs, S_dense) > 0.1
end

@testset "advance! keeps the PDF normalised: $lbl" for (lbl, sde, par, axes, method, mass_tol) in systems
    PI = PathIntegration(sde(copy(par)), method(), Δt, axes()...)
    S = Matrix(PI.stepMX[1])
    q = copy(vec(PI.pdf.p))
    for n in 1:5
        advance!(PI)
        q = S * q
        q ./= PIM._integrate(reshape(q, size(PI.pdf.p)), PI.pdf.axes...)
        @test vec(PI.pdf.p) ≈ q
        @test integrate(PI.pdf) ≈ 1
        @test PI.t ≈ n * Δt
    end
end

@testset "A step approximately conserves probability: $lbl" for (lbl, sde, par, axes, method, mass_tol) in systems
    PI = PathIntegration(sde(copy(par)), method(), Δt, axes()...)
    advance_till_converged!(PI; Tmax = 20.0)
    Sp = reshape(Matrix(PI.stepMX[1]) * vec(PI.pdf.p), size(PI.pdf.p))
    @test abs(1 - PIM._integrate(Sp, PI.pdf.axes...)) < mass_tol # before renormalisation
end

@testset "recompute_PI! with new parameters equals a new PathIntegration: $lbl" for (lbl, sde, par, axes, method, mass_tol) in systems
    par_new = 0.8 .* par
    PI = PathIntegration(sde(copy(par)), method(), Δt, axes()...)
    advance!(PI)
    h(x...) = exp(-sum(abs2, x))
    recompute_PI!(PI; par = par_new, Q_reinit_pdf = true, f = h)
    PI_new = PathIntegration(sde(copy(par_new)), method(), Δt, axes()...)
    @test Matrix(PI.stepMX[1]) ≈ Matrix(PI_new.stepMX[1]) rtol = 1e-10
    @test PI.pdf.p ≈ [h(x...) for x in each_latticecoordinate(PI.pdf)]
    @test PI.t == 0
    @test PI.step_idx == 0

    reinit_PI_pdf!(PI) # default initial PDF
    @test PI.pdf.p ≈ PI_new.pdf.p
end

end
