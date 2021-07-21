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
myp = PDGrid(sde, xs;)

@time tt = PathIntegrationProblem(sde,0.25,xs, precompute=true)#, method = RKMaruyama());

@time begin
    ev = eigs(tt.tpdMX)[2][:,1] .|> real
    myp_ev = deepcopy(myp);
    myp_ev.p .= ev;
    myp_ev.p .= myp_ev.p ./ PathIntegrationMethod._integrate(myp_ev.p, xs)
end

@time for _ in 1:1000
    advance!(tt)
end
res = tt.pdgrid;

begin
    figure(1); clf()
    _x = LinRange(xs[1],xs[end],101)
    plot(_x,myp.(_x),label="Initial PDF" )
    plot(_x,res.(_x),label="Iteration" )
    plot(_x,p_AN(_x), label = "Reference")
    plot(_x,myp_ev.(_x), label = "Eigenvector")
    # scatter(xs,myp_ev.p)
    
    legend()
    # # scatter(xdef,res.(xdef))
end


tt.idx₁ .= [16]
vs = LinRange(-6,6,1001)
_tps = [tt(tt.temp,v) for v in vs]
sum(_tps)*(vs[2]-vs[1])
begin
    figure(1); clf();
    plot(vs,_tps)
end
