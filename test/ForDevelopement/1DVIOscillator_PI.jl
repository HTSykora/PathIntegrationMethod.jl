using Revise, BenchmarkTools
using PathIntegrationMethod

using BenchmarkTools

function f2(x,p,t)
    ζ, _ = p # ζ, σ
    -2ζ*x[2] - x[1]
end
function g2(x,p,t)
    p[2] # = σ
end
par = [0.1,0.1]
##
# method = RK2()

W = Wall(x->0.7-0.05x,-1.)
W = Wall(0.7,-1.)
sde = SDE_VIO(f2,g2,W, par)

axisgrid = (QuinticAxis(-1.,3.,51), QuinticAxis(-3.,3.,51))

# @run sdestep = SDEStep(sde, method, 0.1);
# @run PathIntegrationMethod.compute_velocities_to_impact!(step2, sdestep.ID[])

method = Euler()
Δt = 0.01
##
PI = PathIntegration(sde, method, Δt, axisgrid...);
##
@run PathIntegration(sde, method, Δt, axisgrid...);
