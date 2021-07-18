using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK
##

function fx(u,p,t)
    ζ, _ = p # ζ, σ
    -2ζ*u[2] - u[1]
end
function gx(u,p,t)
    p[2] # = σ
end
sde = SDE_Oscillator1D(fx,gx,par = [0.15, 0.5])
xvs = [Axis(-6,6,31,interpolation=:chebyshev),Axis(-6,6,31,interpolation=:chebyshev)]#,Axis(-6,6,25)]
@time tt = PathIntegrationProblem(sde,0.005,xvs...; precompute=true);

@time for _ in 1:1
    advance!(tt)
end

ev = eigs(tt.tpdMX,nev=2)
# pev1 = ev[2][:,1] .|> real; pev1 ./= sum(pev1)
# pev2 = ev[2][:,2] .|> real; pev2 ./= sum(pev2)
begin
    figure(1); clf()

    # ax = axes(projection="3d")
    res = tt.pdgrid
    X = [res.xs[1][i] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
    Y = [res.xs[2][j] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
    # Data for a three-dimensional line
    scatter3D(X, Y, tt.pdgrid.p)
    X_itp = LinRange(-1.,1.,11)
    Y_itp = [1. for _ in X_itp]
    @time p_itp = tt.pdgrid.(X_itp,Y_itp)
    scatter3D(X_itp, Y_itp, p_itp)
    xlabel(L"x")
    ylabel(L"v")
    zlabel(L"p(x,v)")
    # scatter3D(X, Y, reshape(pev1,size(ttc.pdgrid.p)...))
    # scatter3D(X, Y, reshape(pev2,size(ttc.pdgrid.p)...))
    
    # scatter3D(xvs[2], -xvs[2].*ttc.Δt, zeros(length(xvs[2])))
    # scatter3D([1.1515], [-0.5], ttc.pdgrid(1.1515,-0.5))
end