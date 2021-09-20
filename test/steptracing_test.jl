using Revise
using PathIntegrationMethod
using BenchmarkTools
##

f1(x,p,t) = x[2]^2 + x[1]^2 + x[3]^2 + t
function f2(x,p,t)
    ζ,λ = p
    - x[1] - λ*x[1]^3 - 2*ζ*x[2] + sin(t) + λ*x[3]^3 + 2ζ*x[3]
end
function f3(x,p,t)
    ζ,λ = p
    - x[1] - λ*x[1]^3 - 2*ζ*x[2] + cos(t) - λ*x[3]^2
end
function g(x,p,t)
    sqrt(2)
end

method = Euler()
##
x0 = [1.,2.,3.]; t0 = 0.6; t1 = t0 + 0.1;
par = [0.1,1.1]
sde = SDE((f1,f2,f3), g, par)

@time sdestep = SDEStep(sde, method, x0, similar(x0), t0, t1)

@time PathIntegrationMethod.eval_driftstep!(sdestep)
x0_ref = deepcopy(sdestep.x0)
x1_ref = deepcopy(sdestep.x1)

sdestep.x0[1] = x1_ref[1]
sdestep.x0[2] = x1_ref[2]
sdestep.x1[3] = sdestep.x0[3]

##
@btime PathIntegrationMethod.compute_missing_states_driftstep!($sdestep)
reduce(&, x0_ref .≈ sdestep.x0)
reduce(&, x1_ref .≈ sdestep.x1)