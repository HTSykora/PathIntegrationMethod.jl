# Regression check of the step matrices: stores the matrices of a set of test systems (`save`),
# or compares the current ones to the stored ones (`compare`).
# Run from the package folder with: julia --project benchmark/regression_snapshot.jl save|compare
using PathIntegrationMethod, Serialization
const PIM = PathIntegrationMethod
const FILE = joinpath(@__DIR__, "results", "S_baseline.jls")

f(x,p,t) = x[1] - x[1]^3
g(x,p,t) = sqrt(2)
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]
fp(x,p,t) = x[1] - x[1]^3 + p[1]*sin(2π*t/p[2])
par = [0.5, 0.25, 1.0]
duff() = SDE((f1, f2), g2, copy(par))

systems = [
    "d1 cubic RK2" => () -> PathIntegration(SDE(f,g), RK2(), 0.05, CubicAxis(-3.,3.,41)),
    "d1 cubic RK2 dense" => () -> PathIntegration(SDE(f,g), RK2(), 0.05, CubicAxis(-3.,3.,41); stepMXtype = DenseMX()),
    "d1 linear Euler" => () -> PathIntegration(SDE(f,g), Euler(), 0.05, LinearAxis(-3.,3.,41)),
    "d1 cheb RK2" => () -> PathIntegration(SDE(f,g), RK2(), 0.05, ChebyshevAxis(-3.,3.,31)),
    "d1 trig RK4" => () -> PathIntegration(SDE(f,g), RK4(), 0.05, TrigonometricAxis(-3.,3.,31)),
    "d1 quadgk cubic" => () -> PathIntegration(SDE(f,g), RK4(), 0.05, CubicAxis(-3.,3.,21); discreteintegrator = QuadGKIntegrator()),
    "d1 periodic" => () -> PathIntegration(SDE(fp,g,[1.0,0.5]), Euler(), collect(range(0,0.5,length=6)), CubicAxis(-3.,3.,41)),
    "d2 quintic Euler" => () -> PathIntegration(duff(), Euler(), 0.05, QuinticAxis(-4.,4.,15), QuinticAxis(-4.,4.,15)),
    "d2 quintic RK4 dense" => () -> PathIntegration(duff(), RK4(), 0.05, QuinticAxis(-4.,4.,15), QuinticAxis(-4.,4.,15); stepMXtype = DenseMX()),
    "d2 cheb RK4" => () -> PathIntegration(duff(), RK4(), 0.05, ChebyshevAxis(-4.,4.,11), ChebyshevAxis(-4.,4.,11)),
    "d2 mixed Euler" => () -> PathIntegration(duff(), Euler(), 0.05, CubicAxis(-4.,4.,15), ChebyshevAxis(-4.,4.,11)),
    "d2 linear RK2" => () -> PathIntegration(duff(), RK2(), 0.05, LinearAxis(-4.,4.,15), LinearAxis(-4.,4.,15)),
    # mixed axes: sparse-sparse, dense-dense and dense-sparse
    "d2 cubic×quintic RK4" => () -> PathIntegration(duff(), RK4(), 0.05, CubicAxis(-4.,4.,15), QuinticAxis(-4.,4.,17)),
    "d2 cheb×trig RK4" => () -> PathIntegration(duff(), RK4(), 0.05, ChebyshevAxis(-4.,4.,11), TrigonometricAxis(-4.,4.,13)),
    "d2 cheb×quintic RK4" => () -> PathIntegration(duff(), RK4(), 0.05, ChebyshevAxis(-4.,4.,11), QuinticAxis(-4.,4.,17)),
]

Ss = Dict(lbl => [Matrix(S) for S in mk().stepMX] for (lbl, mk) in systems)

if ARGS[1] == "save"
    serialize(FILE, Ss)
    println("saved ", length(Ss), " systems")
else
    ref = deserialize(FILE)
    for (lbl, _) in systems
        a, b = Ss[lbl], ref[lbl]
        err = maximum(maximum(abs, x - y) / maximum(abs, y) for (x, y) in zip(a, b))
        nzdiff = sum(count(!iszero, x) - count(!iszero, y) for (x, y) in zip(a, b))
        println(rpad(lbl, 24), " max rel err = ", err, "   Δnnz = ", nzdiff, "   bitwise = ", a == b)
    end
end
