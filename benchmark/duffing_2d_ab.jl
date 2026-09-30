# Comparison of two versions of the package (e.g. before and after a change) on the Duffing oscillator (d = 2) with uniform
# and mixed axes: step matrix computation (serial and threaded), `integrate` and the evaluation of the PDF on a grid.
# Run it alternately with the two versions, e.g. with a git worktree of the old commit (and a copy of Manifest.toml in it):
#   julia --project=<worktree> -t auto benchmark/duffing_2d_ab.jl
#   julia --project -t auto benchmark/duffing_2d_ab.jl
using PathIntegrationMethod, Printf

# Duffing oscillator: dx = v dt, dv = (-2ζv + x - λx³) dt + σ dW
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]

# (label, axes): uniform sparse and dense, mixed dense-sparse (both orders), mixed sparse-sparse and dense-dense
cases = (("quintic×quintic", () -> (QuinticAxis(-4.,4.,41), QuinticAxis(-4.,4.,41))),
         ("cheb×cheb", () -> (ChebyshevAxis(-4.,4.,25), ChebyshevAxis(-4.,4.,25))),
         ("quintic×cheb", () -> (QuinticAxis(-4.,4.,41), ChebyshevAxis(-4.,4.,25))),
         ("cheb×quintic", () -> (ChebyshevAxis(-4.,4.,25), QuinticAxis(-4.,4.,41))),
         ("cubic×quintic", () -> (CubicAxis(-4.,4.,41), QuinticAxis(-4.,4.,41))),
         ("cheb×trig", () -> (ChebyshevAxis(-4.,4.,25), TrigonometricAxis(-4.,4.,25))))
besttime(f, n) = minimum(@elapsed(f()) for _ in 1:n)
xs = range(-4., 4., length = 201)
ser = SerialRowComputation(); thr = ThreadedRowComputation()
println("Julia threads: ", Threads.nthreads())
for (lbl, axes) in cases
    PI = PathIntegration(SDE((f1, f2), g2, [0.5, 0.25, 1.0]), RK4(), 0.02, axes()...; rowcomputation = ser)
    recompute_stepMX!(PI; rowcomputation = thr)
    advance!(PI)
    t_s = besttime(() -> recompute_stepMX!(PI; rowcomputation = ser), 15)
    a_s = @allocated recompute_stepMX!(PI; rowcomputation = ser)
    t_t = besttime(() -> recompute_stepMX!(PI; rowcomputation = thr), 30)
    t_i = besttime(() -> (for _ in 1:1000; integrate(PI.pdf); end), 5) / 1000
    PI.(xs, xs')
    t_e = besttime(() -> PI.(xs, xs'), 5)
    a_e = @allocated PI.(xs, xs')
    @printf("%-16s S serial %7.2f ms %6.2f MiB | S %d threads %6.2f ms | integrate %.3f μs | PI.(xs, vs') 201² points %6.2f ms %8.2f MiB\n",
        lbl, 1e3t_s, a_s / 2^20, Threads.nthreads(), 1e3t_t, 1e6t_i, 1e3t_e, a_e / 2^20)
end
