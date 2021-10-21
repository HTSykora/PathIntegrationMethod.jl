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
    x[1] - ε*x[1]^3 + σ*x[2]
end
function f2(x,p,t) 
    _, _, μ = p
    -μ*x[2]
end
function g(x,p,t)
    _, _, μ = p
    sqrt(2μ)
end
sde = SDE((f1,f2),g, [1.,2.0,2.0])

##
euler = Euler()
rk4 = RK4()

##
# # Single test run
Δt = 0.00001
# Δt = 0.000002875
gridaxes = (GridAxis(-6, 6, 101, interpolation = :cubic),
            GridAxis(-4, 4, 51, interpolation = :cubic))
@time PI = PathIntegration(sde, euler, Δt, gridaxes..., pre_compute = true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 31, smart_integration = true,int_limit_thickness_multiplier = 8, sparse_stepMX = true, mPDF_IDs = ((1,),(2,)));


Tmax = 0.6;#1.0
@time for _ in 1:Int((Tmax + sqrt(eps(Tmax))) ÷ Δt)
    advance!(PI)
end

# begin
#     figure(1); clf()
#     # ax = axes(projection="3d")
#     X = [gridaxes[1][i] for i in eachindex(gridaxes[1]), j in eachindex(gridaxes[2])]
#     V = [gridaxes[2][j] for i in eachindex(gridaxes[1]), j in eachindex(gridaxes[2])]
#     # Data for a three-dimensional line
#     scatter3D(X, V, PI.pdf.p)
#     xlabel("x")
#     ylabel("z")
# end

@time update_mPDFs!(PI)

begin
    id = 1
    figure(2); clf()
    # ax = axes(projection="3d")
    _x = LinRange(gridaxes[id][1],gridaxes[id][end],101)
    # Data for a three-dimensional line
    plot(_x, PI.marginal_pdfs[id].(_x))
    plot(_x,PathIntegrationMethod.normal1D.(_x))
    plot(_x,zero(_x))
    xlabel(["x","z"][id])
end
