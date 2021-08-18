## computations for debugging

using Revise, BenchmarkTools
using PathIntegrationMethod

using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using Cubature
##

##
identity(x...) = x
foo(x,y) = sin(0.05x*y^2);
foo(x::Vector{Float64}) = foo(x...)

foo1(x) = foo(x,x-1.)

xs = (Axis(0,6,21; interpolation = :chebyshev),Axis(-6,6,21; interpolation = :chebyshev));

_X = [x for x in xs[1], y in xs[2]]
_Y = [y for x in xs[1], y in xs[2]]
fij = [foo(x,y) for x in xs[1], y in xs[2]]

pgr = PDGrid(fij,xs...)

@btime pgr(2.,4)
@time foo(2.,4)

@time (pgr.xs[1].wts * pgr.xs[2].wts') .* fij |> sum
@time hcubature(foo,[0.,-6.],[6.,6.])

quadgk(foo1,0.,6.)

ik = IntegrationKernel(foo,pgr.xs,length.(pgr.xs))
@btime ik(ik.temp,2.,4.)
ik1(val,x) = ik(val,x,x-1.)

quadgk!(ik1,ik.temp,-2.3,6.)
ik.temp |> sum


pol(x,y) = x^2*y;
pol1(x) = pol(x,x-1)
polij = [pol(x,y) for x in xs[1], y in xs[2]]

f2(x) = foo1(x)*pol1(x)
@time quadgk(f2,-2.3,6.)

viktemp = vec(ik.temp);
vpolij= vec(polij)

@time viktemp' * vpolij
@time ik.temp .* polij |> sum
itg(f) = sum(f(xs[1][i])*ik.temp[i,j] for j in eachindex(xs[2]), i in eachindex(xs[1]))
@time itg(pol1)


##
struct GKTester{fT,fxT,xT,nT}
    f::fT
    fx::fxT
    x::xT
    n::nT
end

function GKTester(f,x_prototype = 0.)
    x = Vector{typeof(x_prototype)}(undef,0);
    n = [0];
    fx = Vector{typeof(f(x_prototype))}(undef,0);
    GKTester(f,fx,x,n)
end
function reinit!(g::GKTester)
    g.n[1] = 0
    resize!(g.x,0)
    resize!(g.fx,0)
    g
end

function (g::GKTester)(x)
    push!(g.x,x)
    g.n[1] += 1
    push!(g.fx,g.f(x))
    return g.fx[end]
end

g = GKTester(foo1)
reinit!(g)
@time quadgk(g,0.,6.)

g.fx


####
##
ξ2x(ξ, a, b) = (ξ + 1)*((b - a)/2) + a # local to global
x2ξ(x, a, b) = 2*(x - a)/(b - a) - 1 # global to local
x2θ(x, a, b) = acos(x2ξ(x, a, b)) # global to angular
function chebygrid(n::Integer)
    cos.(π*(n-1:-1:0)/(n-1))
end
function chebygrid(xa, xb, n::Integer)
    ξ2x.(chebygrid(n), xa, xb)
end
##
d = 0.25
x = 0.251
x = 0.249
x = 0.24863047384206832 -0.24802867532861947
x = 0.24692208514878441
# x = 0.2445369001834514
Δt = 0.01;
α = 0.5
vmin, vmax = -1. ,1.


vᵢ = (x-d)/Δt
vᵢᵢ = (d-x)/(α*Δt)
if vᵢ<vᵢᵢ 
    (vᵢ, vᵢᵢ, vmax)
else
    (vmin, vᵢᵢ, vᵢ)
end


chebygrid(-0.25,0.25,41)