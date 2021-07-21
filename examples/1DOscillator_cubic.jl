using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK, Arpack
##

function fx(u,p,t)
    ζ, λ, _ = p # ζ, λ, σ
    -2ζ*u[2] + u[1] - λ*u[1]^3
end
function gx(u,p,t)
    p[3] # = σ
end
sde = SDE_Oscillator1D(fx,gx,par = [0.15, 0.25, sqrt(0.5)])
# sde = SDE_Oscillator1D(fx,gx,par = [0.2, 0.1, 0.2])
xvs = [Axis(-7,7,31,interpolation = :chebyshev),Axis(-7,7,31,interpolation = :chebyshev)]
Δt = 0.005;
@time tt = PathIntegrationProblem(sde,Δt,xvs..., precompute = true);
tt2 = deepcopy(tt)
# tt = deepcopy(tt2)

@time for _ in 1:1000
    advance!(tt)
end
# @save "./JLD/tt_Duffing_101x101_0.005_500step.jld2" tt

begin
    figure(1); clf()
    
    # ax = axes(projection="3d")
    res = tt.pdgrid
    X = [res.xs[1][i] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
    Y = [res.xs[2][j] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
    # Data for a three-dimensional line
    scatter3D(X, Y, tt.pdgrid.p)
    # scatter3D(xvs[1],-1 .+ 0.25*xvs[1].^3,zero(xvs[1]))
    # xlim(left=-6,right=6)
    # ylim(bottom=-6,top=6)
    zlim(bottom=0,top=0.15)
    # scatter3D(X, Y, tt.pdgrid.p_temp)
    # scatter3D(X, Y, tt.pdgrid.p_temp .- tt.pdgrid.p)
    xlabel(L"x")
    ylabel(L"v")
    zlabel(L"p(x,v)")
    # scatter3D(X, Y, reshape(pev1,size(ttc.pdgrid.p)...))
    # scatter3D(X, Y, reshape(pev2,size(ttc.pdgrid.p)...))
    
    # scatter3D(xvs[2], -xvs[2].*ttc.Δt, zeros(length(xvs[2])))
    # scatter3D([1.1515], [-0.5], ttc.pdgrid(1.1515,-0.5))
end


# Visualisation

xmin, xmax, Nx= -5,5, 1001
ymin, ymax, Ny= -5,5, 1001
_X = [x for x in LinRange(xmin,xmax,Nx), y in LinRange(ymin,ymax,Ny)]
_Y = [y for x in LinRange(xmin,xmax,Nx), y in LinRange(ymin,ymax,Ny)]
_P = [tt.pdgrid(x,y) for x in LinRange(xmin,xmax,Nx), y in LinRange(ymin,ymax,Ny)]

begin
    figure(1); clf()
    
    # ax = axes(projection="3d")
    # Data for a three-dimensional line
    plot_surface(_X, _Y, _P, cmap=PyPlot.cm.jet , antialiased=true)
    # xlim(left=-6,right=6)
    # ylim(bottom=-6,top=6)
    zlim(bottom = 0, top = 0.15)
    xlabel(L"x")
    ylabel(L"v")
    zlabel(L"p(x,v)")
    zticks([0.,0.05,0.1,0.15])
    
    # scatter3D(X, Y, reshape(pev1,size(ttc.pdgrid.p)...))
    # scatter3D(X, Y, reshape(pev2,size(ttc.pdgrid.p)...))
    
    # scatter3D(xvs[2], -xvs[2].*ttc.Δt, zeros(length(xvs[2])))
    # scatter3D([1.1515], [-0.5], ttc.pdgrid(1.1515,-0.5))
end
