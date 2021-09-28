using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

# using QuadGK, Arpack
##

# 1D problem:
f(x,p,t) = x[1]-x[1]^3
g(x,p,t) = sqrt(2)
sde = SDE(f,g)

# Analytic form
_p_AN(x,ε = 1.) = exp(x^2/2 - ε*x^4/4)
function p_AN(xs, ε=1.)
    itg = quadgk(x->_p_AN(x,ε) ,xs[1],xs[end])[1]
    _p_AN.(xs,Ref(ε)) ./ itg
end

function get_PI_err(N, Δt; interpolation = :chebyshev, xmin = -3.0, xmax = 3.0, _testN = 10001, _x = LinRange(xmin, xmax, _testN), Tmax = 10.0, method = Euler())
    gridaxis = GridAxis(_x[1],_x[end],N,interpolation = interpolation)
    PI = PathIntegration(sde, method, Δt, gridaxis, pre_compute = true);
    for _ in 1:Int((Tmax + sqrt(eps(Tmax))) ÷ Δt)
        advance!(PI)
    end
    err = sum(abs, PI.pdf.(_x) .- p_AN(_x)) * ((_x[end] - _x[1])/length(_x));
    PI, err
end

get_div(x) = 10 .^(diff(log10.(x)))
get_errconv(err,Δ) = mean(log10.(get_div(err))./ log10.(get_div(Δ)))

##

@time begin
    Δts = [0.1, 0.01, 0.001, 0.0001, 0.00002875]
    itp_ords = 11:2:51

    err_Δt = Vector{Float64}(undef,0)
    err_N = Vector{Float64}(undef,0)
    _x = LinRange(-3.0, 3.0, 10001)

    for Δt in Δts
        _, err = get_PI_err(itp_ords[end], Δt; method = Euler(), interpolation = :chebyshev, _x =_x, Tmax = 10.0)
        push!(err_Δt, err)
    end

    for N in itp_ords
        _, err = get_PI_err(N, Δts[end-1]; method = Euler(), interpolation = :chebyshev, _testN = 10001, _x =_x, Tmax = 10.0)
        push!(err_N, err)
    end
end
# 1
abs(get_errconv(err_Δt,Δts) -1.) < 0.1



##
1
##
# # Single test run
Δt = 0.0001
Δt = 0.00002875
Tmax = 10.0
gridaxis = GridAxis(-3,3,21,interpolation = :linear)
@time PI = PathIntegration(sde, Euler(), Δt, gridaxis, pre_compute = true);
@time for _ in 1:1#Int((Tmax + sqrt(eps(Tmax))) ÷ Δt)
    advance!(PI)
end
# sum(abs, PI.pdf.(_x) .- p_AN(_x)) * ((_x[end] - _x[1])/length(_x))

begin
    figure(1); clf()
    _x = LinRange(gridaxis[1], gridaxis[end], 10001)
    plot(_x,PI.pdf.(_x),label="Iteration" )
    plot(_x,p_AN(_x), label = "Reference")
    legend()
end

######################
# # Performance tests

# dbg_IK, dbg_kwargs = PathIntegration(sde, Euler(), Δt, gridaxis, pre_compute = true, debug_mode = true );
# @time dbg_stepMX = PathIntegrationMethod.initialize_stepMX(eltype(dbg_IK.pdf.p), dbg_IK.t, length(dbg_IK.pdf))
# @btime PathIntegrationMethod.compute_stepMX($dbg_IK, $dbg_kwargs...)


# dbg_zp = zip(dbg_IK.temp.itpVs,(1,))
# @btime PathIntegrationMethod.reduce_tempprod($dbg_zp...)

# @btime PathIntegrationMethod.get_IK_weights!($dbg_IK, $dbg_IK.kwargs...)
# dbg_intlimits =(first(dbg_IK.int_axes)[1],first(dbg_IK.int_axes)[end]);
# dbg_kwargs = PathIntegrationMethod.cleanup_quadgk_keywords(;dbg_IK.kwargs...);
# @btime quadgk!($dbg_IK, $dbg_IK.temp.itpM, $dbg_intlimits...; $dbg_kwargs...)

# global myint = 0
# function foo(val,x)
#     global myint +=1
#     dbg_IK(val,x)
# end
# quadgk!(foo, dbg_IK.temp.itpM, -3.,3.)

# @btime dbg_IK($dbg_IK.temp.itpM, 0.)  