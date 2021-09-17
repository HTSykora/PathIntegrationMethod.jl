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

xs = [LinRange(start,stop, 2(num - 1)) for (start,stop,num) in grid_dat]
f_interpolated = [pdf(x...) for x in Iterators.product(xs...)]
f_true = [f(x...) for x in Iterators.product(xs...)]

@btime pdf(0.1,1.2,1.1)
@time PathIntegrationMethod.interpolate(pdf.p,pdf.axes,0.1,1.2,1.1; idx_it = pdf.idx_it)
@btime PathIntegrationMethod.interpolate($pdf.p,$pdf.axes,0.1,1.2,1.1; idx_it = $pdf.idx_it)

@time PathIntegrationMethod.basefun_vals_safe!(pdf.axes[1].temp,pdf.axes[1].itp,pdf.axes[1],1.)
@time PathIntegrationMethod.basefun_vals_safe!(pdf.axes[2].temp,pdf.axes[2].itp,pdf.axes[2],1.)
@btime PathIntegrationMethod.basefun_vals_safe!($pdf.axes[3].temp,$pdf.axes[3].itp, $pdf.axes[3],$(-1.5))


@btime PathIntegrationMethod.basefun_vals!($pdf.axes[3].temp,$pdf.axes[3].itp, $pdf.axes[3],$(-1.5))


@btime PathIntegrationMethod.find_idx($pdf.axes[3],$(-1.5))


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