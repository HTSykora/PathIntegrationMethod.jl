using Pkg; Pkg.activate();
using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK, Arpack
##

f1(u,p,t) = u[2]
function f2(u,p,t)
    ζ, σ, _ = p # ζ, σ
    -2ζ*u[2] - u[1] + σ*u[3]
end
function f3(u,p,t)
    _, _, μ = p
    -μ*u[3]
end
function g3(x,p,t)
    _, _, μ = p
    sqrt(2μ)
end
##

par = [0.05, 0.1, 2.]; # ζ, σ = p
sde = SDE((f1,f2,f3),g3,par)
xmin = -4; xmax = 4; xN = 35;
vmin = -4; vmax = 4; vN = 35;
zmin = -4; zmax = 4; zN = 35;
gridaxes = (GridAxis(xmin,xmax,xN,interpolation=:cubic),
        GridAxis(vmin,vmax,vN,interpolation=:cubic),
        GridAxis(zmin,zmax,zN,interpolation=:cubic));

Δt = 0.025
method = Euler(); # method  = RK4()
@time PI = PathIntegration(sde, method, [0.,Δt],gridaxes...; pre_compute=true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 21, smart_integration = true,int_limit_thickness_multiplier = 6, mPDF_IDs = ((1,),(2,),(3,)), σ_init = 1.);

@time recompute_stepMX!(PI, t=0.02)
##
Tmax = 100.;
Tmax = 0.25;
@time advance!(PI)
@time for _ in 1:Int((Tmax + sqrt(eps())) ÷ Δt)
        advance!(PI)
end
# pev1 = ev[2][:,1] .|> real; pev1 ./= sum(pev1)
# pev2 = ev[2][:,2] .|> real; pev2 ./= sum(pev2)

@time update_mPDFs!(PI)
begin
    id = 3
    figure(2); clf()
    # ax = axes(projection="3d")
    _x = LinRange(gridaxes[id][1],gridaxes[id][end],101)
    # Data for a three-dimensional line
    plot(_x, PI.marginal_pdfs[id].(_x))
    # plot(_x,PathIntegrationMethod.normal1D.(_x))
    plot(_x,zero(_x))
    xlabel(["x","v","z"][id])
end
begin
    figure(1); clf()
    for id in 1:3
        # ax = axes(projection="3d")
        _x = LinRange(gridaxes[id][1],gridaxes[id][end],101)
        # Data for a three-dimensional line
        plot(_x, PI.marginal_pdfs[id].(_x), label=[L"x",L"v",L"z"][id])
        # plot(_x,PathIntegrationMethod.normal1D.(_x))
    end
    ylabel(L"p")
    legend()
    plot(_x,zero(_x))
end
reinit_PI_pdf!(PI, PathIntegrationMethod.init_DiagonalNormalPDF(gridaxes...; μ_init = nothing, σ_init = 0.5))