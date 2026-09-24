module PeriodicSDETest
using PathIntegrationMethod, Test

# Scalar system with periodic forcing: dx = (x - x³ + A sin(2πt/T)) dt + √2 dW
f(x,p,t) = x[1] - x[1]^3 + p[1]*sin(2π*t/p[2])
g(x,p,t) = sqrt(2)
par = [1.0, 0.5] # A, T
T = par[2]
ts = collect(range(0, T, length = 6)) # 5 time steps per period
new_PI(ts) = PathIntegration(SDE(f, g, copy(par)), Euler(), ts, CubicAxis(-3., 3., 41))
normalise(q, PI) = q ./ sum(PI.pdf.axes[1].wts .* q)

PI = new_PI(ts)
S = [Matrix(S_n) for S_n in PI.stepMX]

@testset "One step matrix per time interval" begin
    @test length(PI.stepMX) == length(ts) - 1
    for n in eachindex(S)
        @test S[n] ≈ Matrix(new_PI(ts[n:n+1]).stepMX[1]) rtol = 1e-10
    end
    @test maximum(abs, S[1] - S[3]) > 1e-2 * maximum(abs, S[1]) # the forcing is present
end

@testset "Time stepping cycles through the period" begin
    q = copy(PI.pdf.p)
    for n in eachindex(S)
        advance!(PI)
        q = normalise(S[n] * q, PI)
        @test PI.pdf.p ≈ q
        @test PI.t ≈ ts[n+1]
    end
    @test PI.step_idx == length(S)

    advance!(PI) # the next period starts with S₁
    q = normalise(S[1] * q, PI)
    @test PI.pdf.p ≈ q
    @test PI.t ≈ T + ts[2]

    PI_nu = new_PI([0.0, 0.05, 0.15, 0.2]) # nonuniform time steps
    for _ in 1:3
        advance!(PI_nu)
    end
    @test PI_nu.t ≈ 0.2
    advance!(PI_nu)
    @test PI_nu.t ≈ 0.25
end

@testset "Periodic steady state" begin
    PI_ss = new_PI(ts)
    _, ϵ = advance_till_converged!(PI_ss; Tmax = 20.0) # ends at the end of a period
    # the convergence check compares two consecutive time steps instead of states one period apart,
    # so it never converges for time-periodic systems and stops at Tmax
    @test_broken length(ϵ) - 1 < 20.0 / T
    p_start = copy(PI_ss.pdf.p)
    for _ in 1:2
        advance!(PI_ss)
    end
    err_mid = integrate_diff(PI_ss.pdf, p_start)
    for _ in 3:5
        advance!(PI_ss)
    end
    @test integrate_diff(PI_ss.pdf, p_start) < 1e-6 # p(t + T) = p(t)
    @test err_mid > 1e-2 # but p changes within the period
end

@testset "recompute_stepMX! with new time points equals a new PathIntegration" begin
    PI_rc = new_PI(ts)
    for t_new in ([0.0, 0.1, 0.2, 0.3], # fewer time steps
                  collect(0.0:0.1:0.6), # more time steps
                  0.0:0.25:0.5,         # range
                  0.1)                  # single time step
        recompute_stepMX!(PI_rc; t = t_new)
        PI_new = new_PI(t_new)
        @test PI_rc.ts == PI_new.ts
        @test length(PI_rc.stepMX) == length(PI_new.stepMX)
        @test all(Matrix(S1) ≈ Matrix(S2) for (S1, S2) in zip(PI_rc.stepMX, PI_new.stepMX))
    end
end

end
