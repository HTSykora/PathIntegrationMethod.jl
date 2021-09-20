using Revise
using PathIntegrationMethod
using BenchmarkTools

function f(x...) 
    sin(norm(x))
end
get_pdf_vals(pdf, xs) = [pdf(x...) for x in Iterators.product(xs...)]
##
itp = :chebyshev

grid_dat = [(-1,1,15),(0,5,15),(-7,-1,21)]
pdf = ProbabilityDensityFunction([GridAxis(start,stop,num; interpolation = itp) for (start, stop, num) in grid_dat]...; f = f)
xs = [LinRange(start,stop, 2(num - 1)) for (start,stop,num) in grid_dat]


@time f_interpolated1 =get_pdf_vals(pdf, xs)
@time f_true1 = [f(x...) for x in Iterators.product(xs...)]

grid_dat = [(-1,1,16), (0,5,17), (-7,-1,25)]
pdf = ProbabilityDensityFunction([GridAxis(start,stop,num; interpolation = itp) for (start, stop, num) in grid_dat]...; f = f)
xs = [LinRange(start,stop, 2(num - 1)) for (start,stop,num) in grid_dat]
@time f_interpolated2 =get_pdf_vals(pdf, xs)
@time f_true2 = [f(x...) for x in Iterators.product(xs...)]

sum(abs2, f_true2 .- f_interpolated2) < sum(abs2, f_true2 .- f_interpolated2)


# Performance checks
@btime pdf(0.1,1.2,1.1)
@btime PathIntegrationMethod.interpolate($pdf.p,$pdf.axes,0.1,1.2,1.1; idx_it = $pdf.idx_it)

@btime PathIntegrationMethod.basefun_vals_safe!.($pdf.axes,(0.1,1.2,1.1))