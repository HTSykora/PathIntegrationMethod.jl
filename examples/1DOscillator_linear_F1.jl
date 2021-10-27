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
par = [0.05, 0.5, 5.0]# ζ, σ = p
sde = SDE((f1,f2,f3),g3,par)
xmin = -4; xmax = 4; xN = 21;
vmin = -4; vmax = 4; vN = 21;
zmin = -4; zmax = 4; zN = 21;
gridaxes = (GridAxis(xmin,xmax,xN,interpolation=:quintic),
        GridAxis(vmin,vmax,vN,interpolation=:quintic),
        GridAxis(zmin,zmax,zN,interpolation=:quintic));

Δt = 0.01
method = Euler(); # method  = RK4()
@time PI = PathIntegration(sde, method, [0.,Δt],gridaxes...; pre_compute=true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 21, smart_integration = true,int_limit_thickness_multiplier = 6, mPDF_IDs = ((1,),(2,),(3,)), σ_init = 1.);

@time recompute_stepMX!(PI, t=0.02)
##
Tmax = 100.;
Tmax = 25;
@time advance!(PI)
@time for _ in 1:Int((Tmax + sqrt(eps())) ÷ Δt)
        advance!(PI)
end
# pev1 = ev[2][:,1] .|> real; pev1 ./= sum(pev1)
# pev2 = ev[2][:,2] .|> real; pev2 ./= sum(pev2)

@time update_mPDFs!(PI)
# begin
#     id = 3
#     figure(2); clf()
#     # ax = axes(projection="3d")
#     _x = LinRange(gridaxes[id][1],gridaxes[id][end],101)
#     # Data for a three-dimensional line
#     plot(_x, PI.marginal_pdfs[id].(_x))
#     # plot(_x,PathIntegrationMethod.normal1D.(_x))
#     plot(_x,zero(_x))
#     xlabel(["x","v","z"][id])
# end
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
# reinit_PI_pdf!(PI, PathIntegrationMethod.init_DiagonalNormalPDF(gridaxes...; μ_init = nothing, σ_init = 0.5))

##
### Analyitical reference solution
using Pkg; Pkg.activate();
using LinearAlgebra, SparseArrays, StaticArrays
using QuadGK, ZChop, Distributions
using PathIntegrationMethod

using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

function σ2M(p)
    ζ, σ, μ = p;
    σ2 = zeros(3,3); σ2[end,end] = 2μ;
    SMatrix{3,3}(σ2)
end
function AMX(p)
    ζ, σ, μ = p;
    A = SMatrix{3,3}([0. 1. 0.; -1. -2ζ σ; 0. 0. -μ]); 
    #factorize(A)
end

struct DriftGen{AT}
    A::AT
end
function (A::DriftGen)(t)
    exp(A.A*t)
end
DriftGen(p::AbstractVector) = DriftGen(AMX(p))
struct DiffGen{AT,σT,tT}
    A::AT
    σ2::σT
    t::tT
end
function (D::DiffGen)(s)
    exp(D.A*(D.t[1]-s))*D.σ2*exp(D.A'*(D.t[1]-s))
end
function get_Diffval(D::DiffGen,t)
    D.t[1] = t;
    quadgk(D,0.,t)[1]
end
DiffGen(p::AbstractVector) = DiffGen(AMX(p),σ2M(p),[0.])
##

par = [0.05, 0.5, 5.0]
Φ = DriftGen(par)
Φ2 = DiffGen(par)
Φ2.t[1] = 1.
@time diff1 = zchop.(get_Diffval(Φ2,100.0),1e-9)
@time diff2 = zchop.(get_Diffval(Φ2,200.0),1e-9)

dist = MvNormal(zeros(3),Symmetric(diff2))
refpdf0(x,v,z) = pdf(dist,[x,v,z])


xmin = -4; xmax = 4; xN = 21;
vmin = -4; vmax = 4; vN = 21;
zmin = -4; zmax = 4; zN = 21;
gridaxes = (GridAxis(xmin,xmax,xN,interpolation=:chebyshev),
        GridAxis(vmin,vmax,vN,interpolation=:chebyshev),
        GridAxis(zmin,zmax,zN,interpolation=:chebyshev));
ref_pdf = InterpolatedFunction(gridaxes..., f = refpdf0)
mpdfs = PathIntegrationMethod.initialise_mPDF(ref_pdf,((1,),(2,),(3,)))

@time PathIntegrationMethod.update_mPDF!.(mpdfs,Ref(ref_pdf));

@time mpdfs[1](0.)
@time mpdfs[2](0.)
@time mpdfs[3](0.)
begin
    figure(2); clf();
    for id in 1:3
        # ax = axes(projection="3d")
        _x = LinRange(gridaxes[id][1],gridaxes[id][end],1001)
        # Data for a three-dimensional line
        plot(_x, mpdfs[id].(_x), label=[L"x",L"v",L"z"][id])
        # plot(_x,PathIntegrationMethod.normal1D.(_x))
    end
    ylabel(L"p")
    legend()
    # plot(_x,zero(_x))
end

