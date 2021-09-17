using Revise
using PathIntegrationMethod
using BenchmarkTools
##

function f(x...) 
    sin(norm(x))
end

grid_dat = [(-1,1,11),(0,5,11),(-7,-1,21)]
itp = :chebyshev
pdf = ProbabilityDensityFunction([GridAxis(start,stop,num; interpolation = itp) for (start, stop, num) in grid_dat]...; f = f)

xs = [LinRange(start,stop, 10(num - 1)) for (start,stop,num) in grid_dat]
Iterators.product(xs...)


method = Euler()
##
x0 = [1.,2.,3.]; t0 = 0.6; t1 = t0 + 0.1;
par = [0.1,1.1]
sde = SDE((f1,f2,f3), g, par)

@time sdestep = SDEStep(sde,method,x0, similar(x0), t0,t1)

@time eval_driftstep!(sdestep)
x0_ref = deepcopy(sdestep.x0)
x1_ref = deepcopy(sdestep.x1)

sdestep.x0[1] = x1_ref[1]
sdestep.x0[2] = x1_ref[2]
sdestep.x1[3] = sdestep.x0[3]

##
@time PathIntegrationMethod.compute_missing_states_driftstep!(sdestep)
reduce(&, x0_ref .≈ sdestep.x0)
reduce(&, x1_ref .≈ sdestep.x1)