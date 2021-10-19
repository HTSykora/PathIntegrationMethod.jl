using Pkg; Pkg.activate();
using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

# using QuadGK, Arpack
##
# 1D problem:
function f1(x,p,t) 
    ε, σ = p
    x[1]-ε*x[1]^3 +σ*x[2]
end
function f2(x,p,t) 
    _, _, μ = p
    -μ*x[2]
end
function g(x,p,t)
    _, _, μ = p
    sqrt(2μ)
end
sde = SDE((f1,f2),g, [1.,2.,1.])

##
euler = Euler()
rk4 = RK4()

##
# # Single test run
Δt = 0.001
# Δt = 0.000002875
Tmax = 1.0
gridaxes = (GridAxis(-4, 4, 71, interpolation = :chebyshev),
            GridAxis(-4, 4, 35, interpolation = :chebyshev))
@time PI = PathIntegration(sde, rk4, Δt, gridaxes..., pre_compute = true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 31, smart_integration = true,int_limit_thickness_multiplier = 10, sparse_stepMX = true);


@time for _ in 1:Int((Tmax + sqrt(eps(Tmax))) ÷ Δt)
    advance!(PI)
end

begin
    figure(1); clf()
    # ax = axes(projection="3d")
    X = [gridaxes[1][i] for i in eachindex(gridaxes[1]), j in eachindex(gridaxes[2])]
    V = [gridaxes[2][j] for i in eachindex(gridaxes[1]), j in eachindex(gridaxes[2])]
    # Data for a three-dimensional line
    scatter3D(X, V, PI.pdf.p)
    xlabel("x")
    ylabel("z")
end