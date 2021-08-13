using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK, Arpack

struct ExcitationFunction{fT,pT}
    f::Function
    p::pT
end
(ef::ExcitationFunction)(p,t) = ef.f(p,t)
(ef::ExcitationFunction)(t) = ef(ef.p,t)
##

function fx(u,p,t)
    f, g = p # f, g, σ
    f(t) + g
end
function gx(u,p,t)
    p[3] # = σ
end


f_per = ExcitationFunction((p,t)-> cos(π*t + p[1]),[0.]) # p = [φ₀]
β = π/3; r = 0.5;
Mg_F = 0.1245/5 * 9.81 # M*g/|F|
g = Mg_F*sin(β) # M*g/|F|*sin(β)
σ = 0.1
p=[f_per, g, σ]
