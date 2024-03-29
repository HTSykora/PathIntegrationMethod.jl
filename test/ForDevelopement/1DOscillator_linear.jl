using Pkg; Pkg.activate()
using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

##
# steady state second moment: σ²/4ζ
function stM2(p)
    ζ, σ = p
    σ^2 / (4ζ)
end
function p_AN(x,v; σ2 = 1.)
    PathIntegrationMethod.normal1D_σ2(0.,σ2, x)*PathIntegrationMethod.normal1D_σ2(0.,σ2, v)
end
f1(x,p,t) = x[2]
function f2(x,p,t)
    ζ, _ = p # ζ, σ
    -2ζ*x[2] - x[1]
end
function g2(x,p,t)
    p[2] # = σ
end

par = [0.05, 0.1]; # ζ, σ = p
sde = SDE((f1,f2),g2,par)
xmin = -2; xmax = 2; xN = 31;
vmin = -2; vmax = 2; vN = 31;
gridaxes = (ChebyshevAxis(xmin,xmax,xN),
        ChebyshevAxis(vmin,vmax,vN))
# gridaxes = (AxisGrid(xmin,xmax,xN,interpolation=:chebyshev),
#         AxisGrid(vmin,vmax,vN,interpolation=:chebyshev))
Δt = 0.025
method = Euler(); method  = RK4()
# @time PI = PathIntegration(sde, method, [0.,Δt],gridaxes...; pre_compute=true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 21, smart_integration = true,int_limit_thickness_multiplier = 6);
PI = PathIntegration(sde, method, [0.,Δt],gridaxes...; pre_compute=true, discreteintegrator = GaussLegendreIntegrator(), di_N = 21, smart_integration = true,int_limit_thickness_multiplier = 6);

@time recompute_stepMX!(PI,t=0.02)
##
Tmax = 100.;
@time for _ in 1:Int((Tmax + sqrt(eps())) ÷ Δt)
        advance!(PI)
end
# pev1 = ev[2][:,1] .|> real; pev1 ./= sum(pev1)
# pev2 = ev[2][:,2] .|> real; pev2 ./= sum(pev2)
begin
    figure(1); clf()
    σ2 = stM2(par)
    # ax = axes(projection="3d")
    res = PI.pdf
    X = [res.axes[1][i] for i in eachindex(res.axes[1]), j in eachindex(res.axes[2])]
    V = [res.axes[2][j] for i in eachindex(res.axes[1]), j in eachindex(res.axes[2])]
    P_AN = [p_AN(x,v; σ2 = σ2) for x in res.axes[1], v in res.axes[2]]
    # Data for a three-dimensional line
    scatter3D(X, V, res.p)
    scatter3D(X,V,P_AN)

    # X_itp = LinRange(-1.,1.,11)
    # V_itp = [1. for _ in X_itp]
    # @time p_itp = res.(X_itp,Y_itp)
    # scatter3D(X_itp, Y_itp, p_itp)
    # xlabel(L"x")
    # ylabel(L"v")
    # zlabel(L"p(x,v)")
    # scatter3D(X, Y, reshape(pev1,size(ttc.pdgrid.p)...))
    # scatter3D(X, Y, reshape(pev2,size(ttc.pdgrid.p)...))
    
    # scatter3D(xvs[2], -xvs[2].*ttc.Δt, zeros(length(xvs[2])))
    # scatter3D([1.1515], [-0.5], ttc.pdgrid(1.1515,-0.5))
end

