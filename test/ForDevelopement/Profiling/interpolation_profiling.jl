using Pkg; Pkg.activate()
using Revise
using PathIntegrationMethod
using Profile
using BenchmarkTools
using StaticArrays
using PProf

function f(x...) 
    sin(norm(x))
end
get_pdf_vals(f_itp, xs) = [f_itp(x...) for x in Iterators.product(xs...)]

function tester(f,a::Vararg{<:T,N}) where {T,N}
    b = f(a...)
    2b
end

function tester2(f,a::Vararg{Number}) where {T<:Number,N}
    b = f(a...)
    2b
end
function tester3(f,a::Vararg{Number})
    b = f(a...)
    b
end
##

# axisf = LinearAxis
# axisf = CubicAxis
axisf = QuinticAxis
# axisf = ChebyshevAxis
# axisf = TrigonometricAxis


grid_dat = [(QuinticAxis,Float64, -10,10,21),(QuinticAxis,Float64, -10,10,25),(QuinticAxis,Float64,-10,10,31)]

f_itp = InterpolatedFunction(Float64,Tuple(axisf(start,stop,num; xT = et, wT = et) for (axisf,et, start, stop, num) in grid_dat)...; f = f)
##
@time f_itp(1.,2.,3.)

a=[1.,1.,1.];
a=MVector(1.,1,1.);
@btime tester($f_itp,$a...)


Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate = 1 tester(f_itp,a...)
PProf.Allocs.pprof(from_c = false)

Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate = 1 f_itp(1.,1.,1.)
PProf.Allocs.pprof(from_c = false)