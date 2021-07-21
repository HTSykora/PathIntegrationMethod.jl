using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK, Arpack
##

function fx(u,p,t)
    ζ, _ = p # ζ, σ
    -2ζ*u[2] - u[1]
end
function gx(u,p,t)
    p[2] # = σ
end
sde = SDE_Oscillator1D(fx,gx,par = [0.01, 0.1])
##
using PathIntegrationMethod: _tp

_x = [0.,0.]
_x0= [0.,1.]
Δt = 1

xs = LinRange(-6,6,1001)
tox = [_tp(sde,sde.par,_x[1],x,Δt,_x0...,0.) for x in xs]

begin
    figure(1); clf()
    plot(xs,tox)
end