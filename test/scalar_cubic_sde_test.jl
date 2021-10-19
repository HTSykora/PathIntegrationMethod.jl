using Pkg; Pkg.activate();
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

function get_PI_err(N, Δt; interpolation = :chebyshev, xmin = -3.0, xmax = 3.0, _testN = 10001, _x = LinRange(xmin, xmax, _testN), Tmax = 10.0, method = Euler(), kwargs...)
    gridaxis = GridAxis(_x[1],_x[end],N,interpolation = interpolation)
    PI = PathIntegration(sde, method, Δt, gridaxis, pre_compute = true;kwargs...);
    for _ in 1:Int((Tmax + sqrt(eps(Tmax))) ÷ Δt)
        advance!(PI)
    end
    err = sum(abs, PI.pdf.(_x) .- p_AN(_x)) * ((_x[end] - _x[1])/length(_x));
    PI, err
end

get_div(x) = 10 .^(diff(log10.(x)))
get_errconv(err,Δ) = mean(log10.(get_div(err))./ log10.(get_div(Δ)))

##
euler = Euler()
rk4 = RK4()
@time begin
    Δts = [0.1, 0.01, 0.001, 0.0001, 0.00002875]
    itp_ords = 11:2:51

    err_Δt_euler = Vector{Float64}(undef,0)
    err_Δt_rk4 = Vector{Float64}(undef,0)
    err_N = Vector{Float64}(undef,0)
    _x = LinRange(-3.0, 3.0, 10001)

    for Δt in Δts
        _, err = get_PI_err(itp_ords[end], Δt; method = euler, interpolation = :chebyshev, _x =_x, Tmax = 10.0,discreteintegrator = ClenshawCurtisIntegrator(),di_mul = 100)
        push!(err_Δt_euler, err)
    end
    for Δt in Δts
        _, err = get_PI_err(itp_ords[end], Δt; method = rk4, interpolation = :chebyshev, _x =_x, Tmax = 10.0,discreteintegrator = ClenshawCurtisIntegrator(),di_mul = 100)
        push!(err_Δt_rk4, err)
    end

    for N in itp_ords
        _, err = get_PI_err(N, Δts[end-1]; method = method, interpolation = :chebyshev, _testN = 10001, _x =_x, Tmax = 10.0)
        push!(err_N, err)
    end
end
# 1
abs(get_errconv(err_Δt_euler,Δts) -1.) < 0.1
abs(get_errconv(err_Δt_rk4,Δts) -1.) < 0.1

begin
    figure(2); clf();

    legend()
    xscale("log")
    yscale("log")
    xlabel(L"\Delta t")
    ylabel(L"\Delta t")
end
##



##
# # Single test run
Δt = 0.0001
Δt = 0.000002875
Tmax = 10.0
gridaxis = GridAxis(-3, 3, 101, interpolation = :chebyshev)
@time PI = PathIntegration(sde, Euler(), Δt, gridaxis, pre_compute = true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 15, smart_integration = true,int_limit_thickness_multiplier = 6, sparse_stepMX = true);
@time for _ in 1:Int((Tmax + sqrt(eps(Tmax))) ÷ Δt)
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

dbg_IK, dbg_kwargs = PathIntegration(sde, Euler(), Δt, gridaxis, pre_compute = true, debug_mode = true );
@time dbg_stepMX = PathIntegrationMethod.initialize_stepMX(eltype(dbg_IK.pdf.p), dbg_IK.t, length(dbg_IK.pdf))
@btime PathIntegrationMethod.compute_stepMX($dbg_IK, $dbg_kwargs...)


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

### Debugging

# (5,)
# ERROR: LoadError: DomainError with -2.056640625:
# integrand produced Inf in the interval (-2.0625, -2.05078125)

# (2,)
# ERROR: LoadError: DomainError with -0.234375:
# integrand produced Inf in the interval (-0.28125, -0.1875)

ID = 550
gridaxis1 = GridAxis(-3, 3, 1001, interpolation = :linear)
@time dbg_IK1, dbg_kwargs1 = PathIntegration(sde, Euler(), Δt, gridaxis1, pre_compute = true, debug_mode = true );
@time dbg_stepMX1 = PathIntegrationMethod.initialize_stepMX(eltype(dbg_IK1.pdf.p), dbg_IK1.t, length(dbg_IK1.pdf))

@time PathIntegrationMethod.update_IK_state_x1!(dbg_IK1,(ID,))
@time PathIntegrationMethod.update_dyn_state_x1!(dbg_IK1, (ID,))

@btime dbg_IK1($dbg_IK1.temp.itpM,gridaxis1[ID] + 0.001)
dbg_IK1.temp.itpM .* dbg_IK1.pdf.p |> sum
@btime dbg_IK1($dbg_IK1.temp.itpM,-2.7)
@btime dbg_IK1($dbg_IK1.temp.itpM,$(gridaxis[2] + 0.01))
dbg_intlimits =(first(dbg_IK1.int_axes)[1],first(dbg_IK1.int_axes)[end]);
dbg_kwargs = PathIntegrationMethod.cleanup_quadgk_keywords(;dbg_IK1.kwargs...);

gridaxis2 = GridAxis(-3, 3, 1001, interpolation = :chebyshev)
@time dbg_IK2, dbg_kwargs2 = PathIntegration(sde, Euler(), Δt, gridaxis2, pre_compute = true, debug_mode = true );
@time dbg_stepMX2 = PathIntegrationMethod.initialize_stepMX(eltype(dbg_IK2.pdf.p), dbg_IK2.t, length(dbg_IK2.pdf))

dbg_IK2.x1 .= dbg_IK1.x1
dbg_IK2.sdestep.x1 .= dbg_IK1.sdestep.x1 

@btime dbg_IK2($dbg_IK2.temp.itpM,gridaxis1[ID] + 0.001)
dbg_IK2.temp.itpM .* dbg_IK2.pdf.p |> sum

@run dbg_IK2(dbg_IK2.temp.itpM,1.875)
@btime dbg_IK2($dbg_IK2.temp.itpM,-2.7)
@btime dbg_IK2($dbg_IK2.temp.itpM,$(gridaxis[2] + 0.01))
dbg_intlimits =(first(dbg_IK2.int_axes)[1],first(dbg_IK2.int_axes)[end]);
dbg_kwargs = PathIntegrationMethod.cleanup_quadgk_keywords(;dbg_IK2.kwargs...);







global myint = 0
function foo(val,x)
    global myint +=1
    dbg_IK(val,x)
end
quadgk!(foo, dbg_IK.temp.itpM, -3.,3.)
@time quadgk!(dbg_IK, dbg_IK.temp.itpM, dbg_intlimits...; dbg_kwargs..., maxevals=1000)
@btime quadgk!($dbg_IK, $dbg_IK.temp.itpM, $dbg_intlimits...; $dbg_kwargs...)

g_IK, $dbg_IK.temp.itpM, $dbg_intlimits...; $dbg_kwargs...)