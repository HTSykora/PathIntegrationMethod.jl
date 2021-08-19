using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK, Arpack
##

# 1D problem:
f(x,p,t) = x-x^3
g(x,p,t) = sqrt(2)

# Analytic form of the 
_p_AN(x,ε = 1.) = exp(x^2/2 - ε*x^4/4)
function p_AN(xs, ε=1.)
    itg = quadgk(x->_p_AN(x,ε) ,xs[1],xs[end])[1]
    _p_AN.(xs,Ref(ε)) ./ itg
end
sde = SDE(f,g)
xs = Axis(-5,5,31,interpolation = :chebyshev)

@time pip = PathIntegrationProblem(sde,0.25,xs, precompute=true)#, method = RKMaruyama());

@time begin
    ev = eigs(pip.tpdMX)[2][:,1] .|> real
    pip_ev = PDGrid(sde, xs;);
    pip_ev.p .= ev;
    pip_ev.p .= myp_ev.p ./ PathIntegrationMethod._integrate(myp_ev.p, xs)
end

@time for _ in 1:1000
    advance!(pip)
end

begin
    figure(1); clf()
    _x = LinRange(xs[1],xs[end],101)
    plot(_x,myp.(_x),label="Initial PDF" )
    plot(_x,pip.pdgrid.(_x),label="Iteration" )
    plot(_x,p_AN(_x), label = "Reference")
    plot(_x,pip_ev.(_x), label = "Eigenvector")
    
    legend()
end
