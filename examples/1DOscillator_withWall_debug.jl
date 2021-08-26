include("./example_init.jl")
##

function fx(u,p,t)
    ζ = p[1] # f, g, σ
    -2ζ*u[2] - u[1]
end
function gx(u,p,t)
    p[2] # = σ
end
##

ζ = 0.1; σ = 1.
p=[ζ, σ]

r = 0.7
# walls = (Wall(r,-0.1),Wall(r,3.))
walls = (Wall(r,0.),)
sde = SDE_Oscillator1D(fx,gx,par = p)
vi_sde = SDE_VI_Oscillator1D(sde,walls)

Nᵥ = 101; Nₓ = 101; # 300 seconds at this resolution
x_lims = (walls[1].pos, length(walls) >1 ? walls[2].pos : 3.)
x_ax = Axis(x_lims...,Nₓ)
v_ax = Axis(-6,6,Nᵥ)
Δt = 0.05


IK = PathIntegrationProblem(vi_sde, Δt, x_ax, v_ax; precompute=true, μ_init=[2.,0.], σ_init = [0.1,0.25], debug_mode = true);

IK
idx₁ = (5,64)
IK.pdgrid.xs[1][idx₁[1]]
IK.pdgrid.xs[2][idx₁[2]]

PathIntegrationMethod.update_idx1!(IK,idx₁)
PathIntegrationMethod.get_integ_limits(IK)
IK.impactinterval.Q_atwall[1]
IK.impactinterval.lims |> x->(x...,)

@run IK(IK.temp,-2.3)
@run IK(IK.temp,-0.338)


IK.pdgrid.xs[1][IK.idx₁[1]]
IK.pdgrid.xs[2][IK.idx₁[2]]
v₀*Δt₁+ξ
-r(0.)*Δt₂*v₀
get_ξ(IK.method,IK.sde.osc1D, get_t1(IK),get_t0(IK),IK.pdgrid.xs[1][IK.idx₁[1]],nothing,nothing,v₀)
