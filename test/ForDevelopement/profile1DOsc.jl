using Pkg; Pkg.activate()
using Revise, BenchmarkTools
using PathIntegrationMethod
using Profile,ProfileView
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK, Arpack

function fx(u,p,t)
    ζ, λ, _ = p # ζ, λ, σ
    -2ζ*u[2] + u[1] - λ*u[1]^3
end
function gx(u,p,t)
    p[3] # = σ
end
## 
xmin, xmax, Nx = -7, 7, 31
vmin, vmax, Nv = -7, 7, 31

xvs = [Axis(-7,7,Nx,interpolation = :chebyshev),Axis(-7,7,Nv,interpolation = :chebyshev)]

par = [0.15, 0.25, sqrt(0.5)];
# par = [0.2, 0.1, sqrt(0.2)];
# par = [0.2, 0.1, sqrt(0.8)]
# par = [0.02, 0.1, 0.2]
sde = SDE_Oscillator1D(fx,gx,par = par)
Δt = 0.001;
@time pip = PathIntegrationProblem(sde,Δt,xvs..., precompute = true, σ_init = 0.5);


Profile.clear();
@profile PathIntegrationProblem(sde,Δt,xvs..., precompute = true, σ_init = 0.5);
ProfileView.view()

Profile.clear();
##
function tester(a1,a2,n)
    res = PathIntegrationMethod.find_idx(a1,a2[1])
    for i in 1:n
        res = PathIntegrationMethod.find_idx(a1,a2[i])
    end
    res
end
a1 = LinRange(0.,10.,50) |> collect
a2= LinRange(-0.1,10.1,100000)
tester(a1,a2,length(a2))

Profile.clear();
@profile tester(a1,a2,length(a2))
ProfileView.view()
##