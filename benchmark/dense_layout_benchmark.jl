# Layout of dense step matrices: S (gemv 'N' in advance!) or transpose(Sᵀ) (gemv 'T', the rows of S are the contiguous columns of Sᵀ)
# 1. matrix-vector products of random matrices, 2. dense PathIntegrations end-to-end
# Run from the package folder with: julia --project -t auto benchmark/dense_layout_benchmark.jl
using PathIntegrationMethod
besttime(f, n) = minimum(@elapsed(f()) for _ in 1:n)

println("Matrix-vector products of dense matrices (Julia threads: ", Threads.nthreads(), ")")
n_blas = BLAS.get_num_threads()
for nb in (n_blas, 1)
    BLAS.set_num_threads(nb)
    println("BLAS threads $nb")
    for n in (61, 961, 1681, 3721, 6561)
        S = rand(n, n); St = permutedims(S); x = rand(n); y = similar(x); y2 = similar(x)
        mul!(y, S, x); mul!(y2, transpose(St), x)
        reps = max(10, round(Int, 2e8 / n^2))
        tN = besttime(() -> mul!(y, S, x), reps)
        tT = besttime(() -> mul!(y2, transpose(St), x), reps)
        println("  n = $(rpad(n, 5)) S*x (gemv N): $(rpad(round(1e6tN, sigdigits = 3), 7)) μs   transpose(Sᵀ)*x (gemv T): $(rpad(round(1e6tT, sigdigits = 3), 7)) μs   ratio T/N = $(round(tT/tN, digits = 2))   max diff $(maximum(abs, y - y2))")
    end
end
BLAS.set_num_threads(n_blas)

f(x,p,t) = x[1] - x[1]^3
g(x,p,t) = sqrt(2)
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]
fp(x,p,t) = x[1] - x[1]^3 + p[1]*sin(2π*t/p[2])
cases = (("d=1 Chebyshev N=61", () -> PathIntegration(SDE(f,g), RK4(), 0.01, ChebyshevAxis(-3.,3.,61))),
         ("d=2 Chebyshev N=31", () -> PathIntegration(SDE((f1,f2),g2,[0.5,0.25,1.0]), RK4(), 0.01, ChebyshevAxis(-4.,4.,31), ChebyshevAxis(-4.,4.,31))),
         ("d=2 Chebyshev N=61", () -> PathIntegration(SDE((f1,f2),g2,[0.5,0.25,1.0]), RK4(), 0.01, ChebyshevAxis(-4.,4.,61), ChebyshevAxis(-4.,4.,61))),
         ("d=2 quintic N=21 DenseMX", () -> PathIntegration(SDE((f1,f2),g2,[0.5,0.25,1.0]), RK4(), 0.01, QuinticAxis(-4.,4.,21), QuinticAxis(-4.,4.,21); stepMXtype = DenseMX())),
         ("periodic d=1 Chebyshev N=61", () -> PathIntegration(SDE(fp,g,[1.0,0.5]), RK4(), collect(range(0,0.5,length=11)), ChebyshevAxis(-3.,3.,61))))
println("\nDense step matrices: one advance!, advance_till_converged! and steady_state! (Julia threads: ", Threads.nthreads(), ", BLAS threads: ", BLAS.get_num_threads(), ")")
for (lbl, mk) in cases
    PI = mk()
    println(lbl, ": stepMX type ", nameof(typeof(PI.stepMX[1])))
    advance!(PI)
    t_step = besttime(() -> advance!(PI), 200)
    t_conv = besttime(() -> (reinit_PI_pdf!(PI); advance_till_converged!(PI; Tmax = 20.0)), 3)
    line = "   one step $(round(1e6t_step, sigdigits = 3)) μs, advance_till_converged! $(round(1000t_conv, sigdigits = 3)) ms"
    for m in (length(PI.stepMX) == 1 ? (:eigen, :lu, :arnoldi) : (:eigen, :arnoldi))
        t = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            steady_state!(PI; method = m)
            besttime(() -> steady_state!(PI; method = m), 3)
        end
        line *= ", $m $(round(1000t, sigdigits = 3)) ms"
    end
    println(line)
end
