## computations for debugging

using Revise, BenchmarkTools
using PathIntegrationMethod

using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using Cubature
##
function create_PDGrid(fij,xs...)ű
    PDGrid{0,0,eltype(fij),typeof(xs),typeof(fij),Nothing,Nothing,Nothing,Nothing}(xs,fij,nothing,nothing,nothing,nothing)
end


##
foo(x,y) = sin(0.05x*y^2);

xs = (Axis(-6,6,51; interpolation = :chebyshev),Axis(-6,6,51; interpolation = :chebyshev));

_X = [x for x in xs[1], y in xs[2]]
_Y = [y for x in xs[1], y in xs[2]]
fij = [foo(x,y) for x in xs[1], y in xs[2]]

PDGrid{}()