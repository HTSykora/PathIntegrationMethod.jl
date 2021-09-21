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
axis = GridAxis(-3,3,51,interpolation = :chebyshev)

Δt = 0.01
@time pi = PathIntegration(sde, Euler(), Δt, axis, pre_compute = true);

dbg_IK, dbg_kwargs = PathIntegration(sde, Euler(), Δt, axis, pre_compute = true, debug_mode = true );
@time dbg_stepMX = PathIntegrationMethod.initialize_stepMX(eltype(dbg_IK.pdf.p), dbg_IK.t, length(dbg_IK.pdf))
@btime PathIntegrationMethod.compute_stepMX($dbg_IK, $dbg_kwargs...)

@btime PathIntegrationMethod.update()


function fx(dbg_IK)
    PathIntegrationMethod.basefun_vals_safe!.(dbg_IK.temp.itpVs,dbg_IK.pdf.axes,dbg_IK.sdestep.x0 )
    nothing
end
@btime fx($dbg_IK)
@btime PathIntegrationMethod.fill_vals!($dbg_IK.temp.itpM,$dbg_IK,1.,$dbg_IK.temp.idx_it)
dbg_zp = zip(dbg_IK.temp.itpVs,(1,))
@btime PathIntegrationMethod.reduce_tempprod($dbg_zp...)

@btime PathIntegrationMethod.get_IK_weights!($dbg_IK, $dbg_IK.kwargs...)
dbg_intlimits = PathIntegrationMethod.get_int_limits(dbg_IK);
dbg_kwargs = PathIntegrationMethod.cleanup_quadgk_keywords(;dbg_IK.kwargs...);
@btime quadgk!($dbg_IK, $dbg_IK.temp.itpM, $dbg_intlimits...; $dbg_kwargs...)

global myint = 0
function foo(val,x)
    global myint +=1
    dbg_IK(val,x)
end
quadgk!(foo, dbg_IK.temp.itpM, -3.,3.)

@btime dbg_IK($dbg_IK.temp.itpM, 0.)
# @time begin
#     ev = eigs(pip.tpdMX)[2][:,1] .|> real
#     pip_ev = PDGrid(sde, xs;);
#     pip_ev.p .= ev;
#     pip_ev.p .= myp_ev.p ./ PathIntegrationMethod._integrate(myp_ev.p, xs)
# end

@time for _ in 1:1000
    advance!(pip)
end

begin
    figure(1); clf()
    _x = LinRange(xs[1],xs[end],1001)
    plot(_x,pip.pdgrid.(_x),label="Iteration" )
    # plot(_x,p_AN(_x), label = "Reference")
    # plot(_x,pip_ev.(_x), label = "Eigenvector")
    ylim(top=0.8)
    legend()
end
