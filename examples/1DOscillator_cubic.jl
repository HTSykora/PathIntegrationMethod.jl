using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK, Arpack
##
function _p_an(x,v,p)
    ζ, λ, σ = p    
    C0 = 2ζ/(σ^2/2)
    exp(C0*(0.5*x^2  - 0.5*v^2 - 0.25*λ*x^4));
end

function compute_analyitical(xs::Axis,vs::Axis,par)
    pij = [_p_an(x,v,par) for x in xs, v in vs]
    pd = PDGrid(pij,xs,vs)
    C = _integrate(pd); 
    @. pd.p = pd.p / C
    
    pd
end

function fx(u,p,t)
    ζ, λ, _ = p # ζ, λ, σ
    -2ζ*u[2] + u[1] - λ*u[1]^3
end
function gx(u,p,t)
    p[3] # = σ
end
## 
xmin, xmax, Nx = -7, 7, 51
vmin, vmax, Nv = -7, 7, 51

xvs = [Axis(-7,7,Nx,interpolation = :chebyshev),Axis(-7,7,Nv,interpolation = :chebyshev)]

par = [0.15, 0.25, sqrt(0.5)];
# par = [0.2, 0.1, sqrt(0.2)];
# par = [0.2, 0.1, sqrt(0.8)]
# par = [0.02, 0.1, 0.2]
sde = SDE_Oscillator1D(fx,gx,par = par)
Δt = 0.001;
@time PathIntegrationProblem(sde,Δt,xvs..., precompute = true, σ_init = 0.5);
# @btime pip = PathIntegrationProblem($sde,$Δt,$xvs..., precompute = $true, σ_init = $0.5);

@time for _ in 1:1
    advance!(pip)
end

# @time pip_ev, e_val = get_stationary_by_eigenvectors(pip, nev = 20, ev_id = 1)
pd_analytic = compute_analyitical(xvs..., par)

begin
    figure(1); clf()
    
    res = pip.pdgrid
    X = [res.xs[1][i] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
    Y = [res.xs[2][j] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
   
    scatter3D(X, Y, pip.pdgrid.p)
    scatter3D(X, Y, pd_analytic.p)
    # scatter3D(xvs[1],-1 .+ 0.25*xvs[1].^3,zero(xvs[1]))
    # xlim(left=-6,right=6)
    # ylim(bottom=-6,top=6)
    zlim(bottom=0,top=0.15)
    xlabel(L"x")
    ylabel(L"v")
    zlabel(L"p(x,v)")
end


# High-res, interpolated visualisation
_xmin, _xmax, _Nx= -5,5, 101
_vmin, _vmax, _Nv= -5,5, 101
_X = [x for x in LinRange(_xmin,_xmax,_Nx), v in LinRange(_vmin,_vmax,_Nv)]
_V = [v for x in LinRange(_xmin,_xmax,_Nx), v in LinRange(_vmin,_vmax,_Nv)]
_P = [pip.pdgrid(x,y) for x in LinRange(_xmin,_xmax,_Nx), y in LinRange(_vmin,_vmax,_Nv)]
_PAN = [pd_analytic(x,y) for x in LinRange(_xmin,_xmax,_Nx), y in LinRange(_vmin,_vmax,_Nv)]

begin
    figure(1); clf()
    
    # plot_surface(_X, _V, _P, cmap=PyPlot.cm.jet)
    # plot_surface(_X, _V, _PAN, cmap=PyPlot.cm.jet)
    plot_surface(_X, _V, _P.-_PAN, cmap=PyPlot.cm.jet)
    # xlim(left=_xmin,right=_xmax)
    # ylim(bottom=_vmin,top=_vmax)
    zlim(bottom = 0, top = 0.15)
    xlabel(L"x")
    ylabel(L"v")
    zlabel(L"p(x,v)")
    zticks([0.,0.05,0.1,0.15])
    
end
