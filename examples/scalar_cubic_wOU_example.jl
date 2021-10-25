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
    ε, σ =  p
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
sde = SDE((f1,f2),g, [1.,1.0,4.0])

##
euler = Euler()
rk4 = RK4()

##
# # Single test run
Δt = 0.001
# Δt = 0.000002875
@time gridaxes = (GridAxis(-2, 2, 51, interpolation = :cubic),
            GridAxis(-4, 4, 21, interpolation = :cubic))
@time PI = PathIntegration(sde, euler, Δt, gridaxes..., pre_compute = true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 31, smart_integration = true,int_limit_thickness_multiplier = 8, sparse_stepMX = true, mPDF_IDs = ((1,),(2,)), σ_init = 1.);

f_init = deepcopy(PI.pdf)

# @time reinit_PI_pdf!(PI,f_init)
# @time recompute_stepMX!(PI, par = nothing, Q_reinit_pdf = true, f = f_init)

Tmax = 1.0;#1.0
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

##
# Monte-Carlo simulations
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];
using StochasticDiffEq

function MC_f(du,u,p,t)
    ε, σ, μ =  p
    du[1] = u[1] - ε*u[1]^3 + σ*u[2]
    du[2] = -μ*u[2]
end
function MC_g(du,u,p,t)
    _, _, μ = p
    du[2] = sqrt(2μ)
end

tspan = (0.0,10000.0)
Δt = 0.001
u₀ = [0.,0.];

prob = SDEProblem(MC_f,MC_g,u₀,tspan,[1.,1.0,4.0])

@time solu = solve(prob,EM(),dt=Δt)

# begin
#     figure(1); clf();
#     plot(solu.t,getindex.(solu.u,1))
#     # plot(solu.t,getindex.(solu.u,2))
# end

begin
    figure(1); clf();
    hist(getindex.(solu.u,1)[findfirst(x->x>100.,solu.t):end],bins=50,density = true)
    # plot(solu.t,getindex.(solu.u,2))
end
begin
    id = 1
    figure(1); #clf()
    # ax = axes(projection="3d")
    _x = LinRange(gridaxes[id][1],gridaxes[id][end],101)
    # Data for a three-dimensional line
    plot(_x, PI.marginal_pdfs[id].(_x))
    xlabel(["x","z"][id])
end