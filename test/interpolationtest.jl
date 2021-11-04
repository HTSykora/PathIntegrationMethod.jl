using Pkg; Pkg.activate()
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
itp = :cubic
itp = :quintic

grid_dat = [(Float64, -1,1,11),(Float64, 0,5,15),(Float64, -7,-1,21)]
f_itp = InterpolatedFunction(Float64,[GridAxis(start,stop,num; interpolation = itp, xT = et, wT = et) for (et, start, stop, num) in grid_dat]...; f = f)
xs = [LinRange(start,stop, 2(num - 1)) for (et,start,stop,num) in grid_dat]


@time f_interpolated1 =get_pdf_vals(f_itp, xs)
@time f_true1 = [f(x...) for x in Iterators.product(xs...)]

grid_dat = [(Float64, -1,1,16),(Float64, 0,5,21),(Float64, -7,3,27)]
f_itp = InterpolatedFunction(Float64,[GridAxis(et,start,stop,num; interpolation = itp) for (et, start, stop, num) in grid_dat]...; f = f)
xs = [LinRange(start,stop, 2(num - 1)) for (et,start,stop,num) in grid_dat]
@time f_interpolated2 =get_pdf_vals(f_itp, xs)
@time f_true2 = [f(x...) for x in Iterators.product(xs...)]

sum(abs2, f_true1 .- f_interpolated1)/length(f_true1) > sum(abs2, f_true2 .- f_interpolated2)/length(f_true2)


# Performance checks
@btime f_itp(0.1,1.2,1.1)
@btime PathIntegrationMethod.interpolate($f_itp.p,$f_itp.axes,0.1,1.2,1.1;)

@btime PathIntegrationMethod.basefun_vals_safe!.($f_itp.axes,(0.1,1.2,1.1))

##
# Visual checks
using Revise
using PathIntegrationMethod
using BenchmarkTools
using PyPlot
pygui(true)
function f3(x) 
    sin(x)
end
grid_dat3 = (-1,5,21)
grid_dat3 = (0,2π,21)
f_itp3 = InterpolatedFunction(Float64,GridAxis(grid_dat3...; interpolation = :trigonometric, xT = Float64); f = f3)

start,stop,num = grid_dat3
xs = LinRange(start-1.0,stop+1., 10(num-1)+1)
x_ref = LinRange(start,stop, 100(num-1))
@time f_interpolated3 = f_itp3.(xs, allow_extrapolation = false)
begin
    figure(1); clf()
    plot(f_itp3.axes[1],f_itp3.p,color="red","o", label = "Itp points")
    plot(xs,f_interpolated3,"-", markersize=4, label = "Itp full")
    plot(x_ref,f3.(x_ref), label = "Ref")
    # plot(f_itp3.axes[1],f3.(f_itp3.axes[1]), label= "Ref Itp grid")
    legend()
end
