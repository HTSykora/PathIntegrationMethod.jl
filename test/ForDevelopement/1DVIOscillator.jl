using Revise, BenchmarkTools
using PathIntegrationMethod
function f2(x,p,t)
    ζ, _ = p # ζ, σ
    -2ζ*x[2] - x[1]
end
function g2(x,p,t)
    p[2] # = σ
end
par = [0.1,0.1]
##
method = Euler()

W = Wall(x->x,1.)
sde = SDE_VIO(f2,g2,W, par)

@time sdestep = SDEStep(sde, method, 0.1)