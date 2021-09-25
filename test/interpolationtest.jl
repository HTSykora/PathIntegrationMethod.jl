using Revise
using PathIntegrationMethod
using BenchmarkTools

function f(x...) 
    sin(norm(x))
end
get_pdf_vals(f_itp, xs) = [f_itp(x...) for x in Iterators.product(xs...)]
##

itp = :chebyshev
itp = :linear

grid_dat = [(Float64, -1,1,11),(Float64, 0,5,15),(Float64, -7,-1,21)]
f_itp = InterpolatedFunction(Float64,[GridAxis(et,start,stop,num; interpolation = itp) for (et, start, stop, num) in grid_dat]...; f = f)
xs = [LinRange(start,stop, 2(num - 1)) for (et,start,stop,num) in grid_dat]


@time f_interpolated1 =get_pdf_vals(f_itp, xs)
@time f_true1 = [f(x...) for x in Iterators.product(xs...)]

grid_dat = [(Float64, -1,1,16),(Float64, 0,5,21),(Float64, -7,-1,27)]
f_itp = InterpolatedFunction(Float64,[GridAxis(et,start,stop,num; interpolation = itp) for (et, start, stop, num) in grid_dat]...; f = f)
xs = [LinRange(start,stop, 2(num - 1)) for (et,start,stop,num) in grid_dat]
@time f_interpolated2 =get_pdf_vals(f_itp, xs)
@time f_true2 = [f(x...) for x in Iterators.product(xs...)]

sum(abs2, f_true1 .- f_interpolated1)/length(f_true1) > sum(abs2, f_true2 .- f_interpolated2)/length(f_true2)


# Performance checks
@btime f_itp(0.1,1.2,1.1)
@btime PathIntegrationMethod.interpolate($f_itp.p,$f_itp.axes,0.1,1.2,1.1; idx_it = $f_itp.idx_it)

@btime PathIntegrationMethod.basefun_vals_safe!.($f_itp.axes,(0.1,1.2,1.1))

# Visual checks
using PyPlot
pygui(true)
function f3(x) 
    sin(x)
end
grid_dat3 = [(Float64, -1,5,11)]
f_itp3 = InterpolatedFunction(Float64,[GridAxis(et,start,stop,num; interpolation = itp) for (et, start, stop, num) in grid_dat3]...; f = f3)

xs = [LinRange(start,stop, 2(num-1)) for (et,start,stop,num) in grid_dat3]
x_ref = [LinRange(start,stop, 100(num-1)) for (et,start,stop,num) in grid_dat3]
@btime f_interpolated3 = f_itp3.(xs[1])

begin
    figure(1); clf()
    plot(f_itp3.axes[1],f_itp3,color="red","-o")
    plot(xs[1],f_interpolated3,"-o")
    plot(x_ref[1],f3.(x_ref[1]))
    plot(f_itp3.axes[1],f3.(f_itp3.axes[1]))
end