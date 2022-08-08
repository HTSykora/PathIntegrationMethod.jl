using PathIntegrationMethod

function f(x...) 
    exp(-norm(x)/2)
end
get_pdf_vals(f_itp, xs) = [f_itp(x...) for x in Iterators.product(xs...)]



grid_dat = [(Float64, -10,10,21),(Float64, -10,10,25),(Float64,-10,10,31)]
refaxs = Tuple(QuinticAxis(start,stop,5num; xT = et, wT = et) for (et, start, stop, num) in grid_dat)
err_int = InterpolatedFunction(refaxs...)

##

testresults = BitArray(undef,0)
for (axisf,err_ref) in zip((LinearAxis, CubicAxis, QuinticAxis, ChebyshevAxis, TrigonometricAxis), (1.745, 0.309, 0.191, 3.813, 1.177))
    # axisf = LinearAxis
    f_itp = InterpolatedFunction(Float64,[axisf(start,stop,num; xT = et, wT = et) for (et, start, stop, num) in grid_dat]...; f = f)
    for (i,x) in enumerate(Iterators.product(err_int.axes...))
        err_int.p[i] = abs(f(x...) - f_itp(x...))
    end
    push!(testresults, integrate(err_int) < err_ref)
end

reduce(&,testresults)
