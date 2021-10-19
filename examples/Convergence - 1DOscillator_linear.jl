using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

# using QuadGK, Arpack
##

# 1D problem:
f1(x,p,t) = x[2]
function f2(x,p,t)
    ζ, _ = p # ζ, σ
    -2ζ*x[2] - x[1]
end
function g2(x,p,t)
    p[2] # = σ
end

par = [0.05, 0.1]; # ζ, σ = p
sde = SDE((f1,f2),g2,par)

# Analytic form
function stM2(p)
    ζ, σ = p
    σ^2 / (4ζ)
end
function p_AN(x,v; σ2 = 1.)
    PathIntegrationMethod.normal1D_σ2(0.,σ2, x)*PathIntegrationMethod.normal1D_σ2(0.,σ2, v)
end
function get_PI_err_Δt!(PI,Δt, errF; Tmax = 100., Q_reinit = false) 
    recompute_stepMX!(PI, t=Δt, Q_reinit = Q_reinit)
    for _ in 1:Int((Tmax + sqrt(eps(Tmax))) ÷ Δt)
        advance!(PI)
    end
    
    recycle_interpolatedfunction!(errF, (x,v) -> p_AN(x,v,σ2 = stM2(PI.IK.sdestep.sde.par)))
    errF.p .= abs.(errF.p .- PI.pdf.p)
    PI, integrate(errF)

end


get_div(x) = 10 .^(diff(log10.(x)))
get_errconv(err,Δ) = mean(log10.(get_div(err))./ log10.(get_div(Δ)))

##

xmin = -2; xmax = 2; xN = 51;
vmin = -2; vmax = 2; vN = 51;
gridaxes = (GridAxis(xmin,xmax,xN,interpolation=:chebyshev),
        GridAxis(vmin,vmax,vN,interpolation=:chebyshev))
Δt = 0.025
euler = Euler();
rk4  = RK4();
@time PI_euler = PathIntegration(sde, euler, [0.,Δt],gridaxes...; pre_compute=true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 21, smart_integration = true,int_limit_thickness_multiplier = 6);
@time PI_rk4 = PathIntegration(sde, rk4, [0.,Δt],gridaxes...; pre_compute=true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 21, smart_integration = true,int_limit_thickness_multiplier = 6);

errF = InterpolatedFunction(gridaxes...);

@time recompute_stepMX!(PI_euler,t=0.02)
@time recompute_stepMX!(PI_rk4,t=0.02)

##
@time begin
    Δts = [4.,2.,1., 1/4, 1/64, 1/256, 1/1024,1/4096,1/8192]#0.00002875]

    err_Δt_euler = Vector{Float64}(undef,0)
    err_Δt_rk4 = Vector{Float64}(undef,0)

    for Δt in Δts
        _, err = get_PI_err_Δt!(PI_euler, Δt, errF; Tmax = 100., Q_reinit = false) 
        push!(err_Δt_euler, err)
    end
    for Δt in Δts
        _, err = get_PI_err_Δt!(PI_rk4, Δt, errF; Tmax = 100., Q_reinit = false)
        push!(err_Δt_rk4, err)
    end

    # for N in itp_ords
    #     _, err = get_PI_err(N, Δts[end-1]; method = method, interpolation = :chebyshev, _testN = 10001, _x =_x, Tmax = 10.0)
    #     push!(err_N, err)
    # end
end
# 1
abs(get_errconv(err_Δt_euler,Δts) -1.) < 0.1
abs(get_errconv(err_Δt_rk4,Δts) -1.) < 0.1


begin
    figure(1); clf()
    plot(Δts, err_Δt_euler, label = "Euler, \$N_x = $(xN), N_v = $(vN)\$")
    plot(Δts, err_Δt_rk4, label = "RK4, \$N_x = $(xN), N_v = $(vN)\$")

    plot(Δts, Δts,"-",c = py_colors[8],alpha = 0.5, label=L"\Delta t^{-1}")
    xscale(:log); yscale(:log)
    legend()

    xlabel("\$\\Delta t\$ - Time step")
    ylabel(L"\displaystyle \int \left| p_{\mathrm{ref}}(x) - \tilde{p}_{\Delta t}(x) \right| \mathrm{d}x")
end



@time begin
    Δt = 0.001
    itp_ords = 11:2:51
    itp_ords_pol = 11:20:211

    err_N_cheb = Vector{Float64}(undef,0)
    err_N_lin = Vector{Float64}(undef,0)
    err_N_cub = Vector{Float64}(undef,0)
    _x = LinRange(-3.0, 3.0, 10001)

    for N in itp_ords
        _, err = get_PI_err(N, Δt; method = euler, interpolation = :chebyshev, _testN = 10001, _x =_x, Tmax = 10.0,discreteintegrator = ClenshawCurtisIntegrator(),di_mul = 100)
        push!(err_N_cheb, err)
    end
    for N in itp_ords_pol
        _, err = get_PI_err(N, Δt; method = euler, interpolation = :linear, _testN = 10001, _x =_x, Tmax = 10.0,discreteintegrator = ClenshawCurtisIntegrator(),di_mul = 100)
        push!(err_N_lin, err)
    end
    for N in itp_ords_pol
        _, err = get_PI_err(N, Δt; method = euler, interpolation = :cubic, _testN = 10001, _x =_x, Tmax = 10.0,discreteintegrator = ClenshawCurtisIntegrator(),di_mul = 100)
        push!(err_N_cub, err)
    end
end
begin
    figure(2); clf()
    plot(itp_ords, err_N_cheb,".-",label="Chebyshev \$\\Delta t = $(Δt)\$")
    
    plot(itp_ords_pol, err_N_lin,".-",label="Linear, \$\\Delta t = $(Δt)\$")
    plot(itp_ords_pol, err_N_cub,".-",label="Cubic, \$\\Delta t = $(Δt)\$")
    plot(itp_ords[1:10], exp.(-0.5itp_ords[1:10] .+ 5) ,"-",c = py_colors[8],alpha = 0.5, label = L"\exp(-a N + b)")
    plot(itp_ords_pol, 10itp_ords_pol.^-1 ,"-",c = py_colors[8],alpha = 0.6, label = L"N^{-1}")
    plot(itp_ords_pol, 300itp_ords_pol.^-2 ,"-",c = py_colors[8],alpha = 0.7, label = L"N^{-2}")
    plot(itp_ords_pol, 500itp_ords_pol.^-3 ,"-",c = py_colors[8],alpha = 0.8, label = L"N^{-3}")

    xscale(:log); yscale(:log)

    legend()
    xlabel("\$N\$ - Interpolation order")
    ylabel(L"\displaystyle \int \left| p_{\mathrm{ref}}(x) - \tilde{p}_{N}(x) \right| \mathrm{d}x")
end
# 1




Δt = 0.001
Tmax = 10.
# # Single test run
gridaxis = GridAxis(-3,3,201,interpolation = :linear,newton_cotes_order = 3)
@time PI = PathIntegration(sde, Euler(), Δt, gridaxis, pre_compute = true,discreteintegrator = ClenshawCurtisIntegrator(),di_mul = 100);
@time for _ in 1:1:Int((Tmax + sqrt(eps(Tmax))) ÷ Δt)
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

# @btime dbg_IK($dbg_IK.temp.itpM, 0.)r