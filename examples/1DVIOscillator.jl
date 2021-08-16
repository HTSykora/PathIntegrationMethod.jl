using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK, Arpack

struct ExcitationFunction{fT,pT}
    f::fT
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
##

β = π/18; r = 0.3; d = 0.25;
Mg_F = 0.1245/5 * 9.81 # M*g/|F|
g = Mg_F*sin(β) # M*g/|F|*sin(β)
σ = 0.1
p=[f_per, g, σ]

r = 0.3

sde = SDE_Oscillator1D(fx,gx,par = p)

Nᵥ = 31; Nₓ = 31;
v_ax = Axis(-1,1,31)
ts = LinRange(0.,2.,31)
vi_sde, pdgrid = create_symmetric_VI_PDGrid(sde, d, r, v_ax, Nₓ)

pip = PathIntegrationProblem(vi_sde, pdgrid, ts; precompute=true);

@run PathIntegrationProblem(vi_sde, pdgrid, ts; precompute=true);


IK = PathIntegrationMethod.computeintegrationmatrix(vi_sde, pdgrid, ts, EulerMaruyama())
IK.idx₁ .= [2,1]
@run ik = IK(IK.temp,0.052264231633826985)
quadgk!(IK,IK.temp,v_ax[16],v_ax[17])

PathIntegrationMethod.get_ξ(EulerMaruyama(),vi_sde,ts[2],ts[1],pdgrid.xs[1][2],nothing,nothing,-1)

