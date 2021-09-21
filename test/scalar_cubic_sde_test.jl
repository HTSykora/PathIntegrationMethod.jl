using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

# using QuadGK, Arpack
##

# 1D problem:
f(x,p,t) = x[1]-x[1]^3
g(x,p,t) = sqrt(2)

# Analytic form of the 
_p_AN(x,ε = 1.) = exp(x^2/2 - ε*x^4/4)
function p_AN(xs, ε=1.)
    itg = quadgk(x->_p_AN(x,ε) ,xs[1],xs[end])[1]
    _p_AN.(xs,Ref(ε)) ./ itg
end
sde = SDE(f,g)
gridaxis = GridAxis(-3,3,51,interpolation = :chebyshev)

Δt = 0.01
@time pi = PathIntegration(sde, Euler(), Δt, gridaxis, pre_compute = true);

@time for _ in 1:1000
    advance!(pi)
end

begin
    figure(1); clf()
    _x = LinRange(axis[1],axis[end],1001)
    plot(_x,pi.pdf.(_x),label="Iteration" )
    plot(_x,p_AN(_x), label = "Reference")
    # plot(_x,pip_ev.(_x), label = "Eigenvector")
    # ylim(top=0.8)
    legend()
end


######################
# Performance tests

dbg_IK, dbg_kwargs = PathIntegration(sde, Euler(), Δt, axis, pre_compute = true, debug_mode = true );
@time dbg_stepMX = PathIntegrationMethod.initialize_stepMX(eltype(dbg_IK.pdf.p), dbg_IK.t, length(dbg_IK.pdf))
@btime PathIntegrationMethod.compute_stepMX($dbg_IK, $dbg_kwargs...)


dbg_zp = zip(dbg_IK.temp.itpVs,(1,))
@btime PathIntegrationMethod.reduce_tempprod($dbg_zp...)

@btime PathIntegrationMethod.get_IK_weights!($dbg_IK, $dbg_IK.kwargs...)
dbg_intlimits =(first(dbg_IK.int_axes)[1],first(dbg_IK.int_axes)[end]);
dbg_kwargs = PathIntegrationMethod.cleanup_quadgk_keywords(;dbg_IK.kwargs...);
@btime quadgk!($dbg_IK, $dbg_IK.temp.itpM, $dbg_intlimits...; $dbg_kwargs...)

global myint = 0
function foo(val,x)
    global myint +=1
    dbg_IK(val,x)
end
quadgk!(foo, dbg_IK.temp.itpM, -3.,3.)

@btime dbg_IK($dbg_IK.temp.itpM, 0.)