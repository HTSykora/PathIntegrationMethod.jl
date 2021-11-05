using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];
using JLD2
##

function run_IK_evaluations(IK,idxs)
    for idx in idxs
        PathIntegrationMethod.update_IK_state_x1!(IK, idx)
        PathIntegrationMethod.update_dyn_state_x1!(IK, idx)
        PathIntegrationMethod.rescale_discreteintegrator!(IK; IK.kwargs...)
        PathIntegrationMethod.get_IK_weights!(IK)
    end
    nothing
end

function test_IK(IK, n0 = 1, n1 = 10; ns_scale = 10^6,reducer = median)
    idxs = collect(PathIntegrationMethod.dense_idx_it(IK))[n0:n1];
    bm = @benchmark run_IK_evaluations($IK,$idxs)
    reducer(bm).time /(n1 - n0 + 1)/ns_scale, bm
    # median time in [ms], bemnchmarkresults
end
##
# 1D problem
f(x,p,t) = x[1]-x[1]^3
g(x,p,t) = sqrt(2)
sde = SDE(f,g)
Δt = 0.001

gridaxis = GridAxis(-3, 3, 11, interpolation = :chebyshev)
@time IK = PathIntegration(sde, Euler(), Δt, gridaxis, pre_compute = true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 21, smart_integration = true,int_limit_thickness_multiplier = 6, sparse_stepMX = true, extract_IK = Val{true}());
t,_ = test_IK(IK,1,11)

Ns = 11:10:151
itps = [:chebyshev, :trigonometric, :linear, :cubic, :quintic]
median_times = Vector{Vector{Float64}}(undef,0)
for itp in itps 
    _times = Vector{Float64}(undef,0);
    for N in Ns
        gridaxis = GridAxis(-3, 3, N, interpolation = itp)
        IK = PathIntegration(sde, Euler(), Δt, gridaxis, pre_compute = true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 21, smart_integration = true,int_limit_thickness_multiplier = 6, sparse_stepMX = true, extract_IK = Val{true}());
        t, _ = test_IK(IK,1,11)
        push!(_times, t)
    end
    push!(median_times, _times)
end
@save "./examples/Benchmarks/1_scalar_mediantimes.jld2" median_times itps Ns

begin
    figure(1); clf()
    for (i,_times) in enumerate(median_times)
        plot(Ns,_times,label=itps[i])
    end
    xscale("log"); yscale("log")
    legend()
end
##

# 2D problem
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
xmin = -2; xmax = 2; xN = 21;
vmin = -2; vmax = 2; vN = 21;
Δt = 0.025
gridaxes = (GridAxis(xmin,xmax,xN,interpolation=:chebyshev),
        GridAxis(vmin,vmax,vN,interpolation=:chebyshev))

@time IK = PathIntegration(sde, Euler(), Δt, gridaxes..., pre_compute = true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 21, smart_integration = true,int_limit_thickness_multiplier = 6, sparse_stepMX = true, extract_IK = Val{true}());
t,bm = test_IK(IK,1,11,reducer = mean)

@time PI = PathIntegration(sde, Euler(), Δt, gridaxes..., pre_compute = true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 21, smart_integration = true,int_limit_thickness_multiplier = 6, sparse_stepMX = true);

Ns = 11:10:151
itps = [:chebyshev, :trigonometric, :linear, :cubic, :quintic]
median_times = Vector{Vector{Float64}}(undef,0)
for itp in itps 
    _times = Vector{Float64}(undef,0);
    for N in Ns
        gridaxes = (GridAxis(xmin,xmax,N,interpolation=itp),
                    GridAxis(vmin,vmax,N,interpolation=itp))
        IK = PathIntegration(sde, Euler(), Δt, gridaxes..., pre_compute = true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 21, smart_integration = true,int_limit_thickness_multiplier = 6, sparse_stepMX = true, extract_IK = Val{true}());
        t, _ = test_IK(IK,1,11)
        push!(_times, t)
    end
    push!(median_times, _times)
end
@save "./examples/Benchmarks/2_1Dlinosc_mediantimes.jld2" median_times itps Ns

begin
    figure(2); clf()
    for (i,_times) in enumerate(median_times)
        plot(Ns,_times,label=itps[i])
    end
    xscale("log"); yscale("log")
    legend()
end

##
# 3D problem
f1(u,p,t) = u[2]
function f2(u,p,t)
    ζ, σ, _ = p # ζ, σ
    -2ζ*u[2] - u[1] + σ*u[3]
end
function f3(u,p,t)
    _, _, μ = p
    -μ*u[3]
end
function g3(x,p,t)
    _, _, μ = p
    sqrt(2μ)
end
##
par = [0.05, 0.5, 5.0]# ζ, σ = p
sde = SDE((f1,f2,f3),g3,par)
xmin = -4; xmax = 4; xN = 31;
vmin = -4; vmax = 4; vN = 31;
zmin = -4; zmax = 4; zN = 31;
Δt = 0.025
gridaxes = (GridAxis(xmin,xmax,xN,interpolation=:chebyshev),
        GridAxis(vmin,vmax,vN,interpolation=:chebyshev),
        GridAxis(vmin,vmax,vN,interpolation=:chebyshev))

@time IK = PathIntegration(sde, Euler(), Δt, gridaxes..., pre_compute = true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 21, smart_integration = true,int_limit_thickness_multiplier = 6, sparse_stepMX = true, extract_IK = Val{true}());

Ns = 11:10:151
itps = [:chebyshev, :trigonometric, :linear, :cubic, :quintic]
median_times = Vector{Vector{Float64}}(undef,0)
for itp in itps 
    _times = Vector{Float64}(undef,0);
    for N in Ns
        gridaxes = (GridAxis(xmin,xmax,N,interpolation=itp),
                    GridAxis(vmin,vmax,N,interpolation=itp),
                    GridAxis(vmin,vmax,N,interpolation=itp))
        IK = PathIntegration(sde, Euler(), Δt, gridaxes..., pre_compute = true, discreteintegrator = ClenshawCurtisIntegrator(), di_N = 21, smart_integration = true,int_limit_thickness_multiplier = 6, sparse_stepMX = true, extract_IK = Val{true}());
        t, _ = test_IK(IK,1,11)
        push!(_times, t)
    end
    push!(median_times, _times)
end
@save "./examples/Benchmarks/3_1DlinoscF1_mediantimes.jld2" median_times itps Ns

begin
    figure(3); clf()
    for (i,_times) in enumerate(median_times)
        plot(Ns,_times,label=itps[i])
    end
    xscale("log"); yscale("log")
    legend()
end

##
# Read results
_i = 1;
files = ["./examples/Benchmarks/1_scalar_mediantimes.jld2",
    "./examples/Benchmarks/2_1Dlinosc_mediantimes.jld2",
    "./examples/Benchmarks/3_1DlinoscF1_mediantimes.jld2"]
@load files[_i] 

begin
    figure(_i); clf()
    for (i,_times) in enumerate(median_times)
        plot(Ns,_times,label=itps[i])
    end
    xscale("log"); yscale("log")
    legend()
end