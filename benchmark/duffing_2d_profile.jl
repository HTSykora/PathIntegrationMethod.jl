# Duffing oscillator (d = 2) with uniform (sparse or dense) and mixed dense-sparse axes: timings and allocations of the
# step matrix computation and of the functions used in a time step (`bench`), and a snapshot of the results to check
# that a change does not change them (`save`, `compare`).
# Run from the package folder with: julia --project -t auto benchmark/duffing_2d_profile.jl bench|save|compare
using PathIntegrationMethod, Serialization, Printf
const PIM = PathIntegrationMethod
const FILE = joinpath(@__DIR__, "results", "duffing_2d_snapshot.jls")

# Duffing oscillator: dx = v dt, dv = (-2ζv + x - λx³) dt + σ dW
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]
par = [0.5, 0.25, 1.0] # ζ, λ, σ

# (label, axes)
cases = (("quintic × quintic", () -> (QuinticAxis(-4., 4., 41), QuinticAxis(-4., 4., 41))),
         ("Chebyshev × Chebyshev", () -> (ChebyshevAxis(-4., 4., 25), ChebyshevAxis(-4., 4., 25))),
         ("quintic × Chebyshev", () -> (QuinticAxis(-4., 4., 41), ChebyshevAxis(-4., 4., 25))),
         ("Chebyshev × quintic", () -> (ChebyshevAxis(-4., 4., 25), QuinticAxis(-4., 4., 41))))
serial = SerialRowComputation()
new_PI(axes; kwargs...) = PathIntegration(SDE((f1, f2), g2, copy(par)), RK4(), 0.02, axes()...; mPDF_IDs = [1, 2], kwargs...)

# Minimum time of `n` evaluations of `f()`
besttime(f, n = 5) = minimum(@elapsed(f()) for _ in 1:n)
# time and allocations per call of `f()`, repeated `n` times inside a function
function percall(f, n)
    f()
    t = @elapsed for _ in 1:n; f(); end
    a = @allocated for _ in 1:n; f(); end
    t / n, a / n
end
xvs = [(x, v) for x in range(-3.9, 3.9, length = 21) for v in range(-3.9, 3.9, length = 21)]
eval_grid(PI) = sum(PI(x, v) for (x, v) in xvs)

function results(axes)
    PI = new_PI(axes; rowcomputation = serial)
    S = [Matrix(S) for S in PI.stepMX]
    for _ in 1:20; advance!(PI); end
    update_mPDFs!(PI)
    (; S, wts = copy(PI.stepMX_wts), p = copy(PI.pdf.p), mass = integrate(PI.pdf), Ex2 = integrate((x, v) -> x^2, PI.pdf),
       diff = integrate_diff(PI.pdf, PI.p_temp), vals = [PI(x, v) for (x, v) in xvs], m1 = copy(PI.marginal_pdfs[1].pdf.p), m2 = copy(PI.marginal_pdfs[2].pdf.p))
end

mode = isempty(ARGS) ? "bench" : ARGS[1]
if mode == "bench"
    println("Julia threads: ", Threads.nthreads(), ", BLAS threads: ", PIM.BLAS.get_num_threads())
    for (lbl, axes) in cases
        PI = new_PI(axes; rowcomputation = serial)
        recompute_stepMX!(PI; rowcomputation = ThreadedRowComputation())
        t_S = besttime(() -> recompute_stepMX!(PI; rowcomputation = serial))
        a_S = @allocated recompute_stepMX!(PI; rowcomputation = serial)
        t_Sp = besttime(() -> recompute_stepMX!(PI; rowcomputation = ThreadedRowComputation()))
        t_adv, a_adv = percall(() -> advance!(PI), 200)
        t_int, a_int = percall(() -> integrate(PI.pdf), 200)
        t_diff, a_diff = percall(() -> integrate_diff(PI.pdf, PI.p_temp), 200)
        t_eval, a_eval = percall(() -> eval_grid(PI), 20)
        t_m, a_m = percall(() -> update_mPDFs!(PI), 200)
        @printf("%-22s %4d rows | S serial %.4f s (%.2f MiB), S %d threads %.4f s | advance! %.1f μs (%.0f B) | integrate %.2f μs (%.0f B) | integrate_diff %.2f μs (%.0f B) | PI(x,v) %.0f ns (%.1f B) | update_mPDFs! %.2f μs (%.0f B)\n",
            lbl, length(PI.pdf), t_S, a_S / 2^20, Threads.nthreads(), t_Sp, 1e6t_adv, a_adv, 1e6t_int, a_int, 1e6t_diff, a_diff,
            1e9t_eval / length(xvs), a_eval / length(xvs), 1e6t_m, a_m)
    end
elseif mode == "save"
    serialize(FILE, Dict(lbl => results(axes) for (lbl, axes) in cases))
    println("saved ", length(cases), " cases to ", FILE)
elseif mode == "compare"
    ref = deserialize(FILE)
    for (lbl, axes) in cases
        r, r0 = results(axes), ref[lbl]
        println(rpad(lbl, 22), " bitwise equal: ", join(("$k: $(getfield(r, k) == getfield(r0, k))" for k in keys(r)), ", "))
        all(getfield(r, k) == getfield(r0, k) for k in keys(r)) ||
            println("    max relative differences: ", join(("$k: $(maximum(abs, getfield(r, k) - getfield(r0, k)) / maximum(abs, getfield(r0, k)))" for k in (:S, :p, :vals)), ", "))
    end
end
