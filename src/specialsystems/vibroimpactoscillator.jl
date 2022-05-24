
###################
(w::Wall{<:Function})(v) = w.r(v)
(w::Wall{<:Number})(v) = w.r
Wall(r) = Wall(r,0.);

# sde/sde.jl
get_dkm(sde::SDE_VIO) = (2,2,1) # (d,k,m)
function osc_f1(u,p,t)
    u[2]
end
function SDE_VIO(f::fT,g::gT, wall::Union{Tuple{wT1},Tuple{wT1,wT2}}, par=nothing) where {fT<:Function,gT<:Function, wT1<:Wall, wT2<:Wall}
    sde = SDE((osc_f1,f), g, par);
    SDE_VIO(sde, wall)
end

SDE_VIO(f::fT,g::gT, w::Wall, par = nothing; kwargs...) where {fT<:Function,gT<:Function} = SDE_VIO(f, g, (w,), par; kwargs...)

# sde/sdestep.jl
function NonSmoothSDEStep(sdesteps::Vararg{SDEStep, N}; kwargs...) where {N}
    @assert reduce(&, Q_compatible(sdesteps[1], sdesteps[i]) for i in 2:N) "Incompatible SDE steps provided"
    d,k,m = get_dkm(sdesteps[1])
    NonSmoothSDEStep{d,k,m,typeof(sdesteps)}(sdesteps)

end

function SDEStep(sde::sdeT, method::methodT, x0,x1, t0, t1; precomputelevel::pclT = PreComputeNewtonStep(), kwargs...) where {sdeT<:SDE_VIO, methodT <: DiscreteTimeSteppingMethod, pclT <: PreComputeLevel} where {d,k,m}
    
    _method = DiscreteTimeStepping(sde, method)
    steptracers = precomputelevel(sde,_method, x0,x1, t0, t1)
    ti = Ref(zero(typeof(t0)));
    step1 = SDEStep{2,2,1,sdeT,typeof(_method),typeof(steptracers[1]),typeof(x0),typeof(x1),typeof(t0), Nothing, Nothing}(sde, _method, similar(x0), similar(x1), t0, t1, steptracers[1], nothing, nothing)
    step2 = SDEStep{2,2,1,sdeT,typeof(_method),typeof(steptracers[2]),typeof(x0),typeof(x1),typeof(t0), typeof(ti), typeof(x0)}(sde, _method, similar(x0), similar(x1), t0, t1, steptracers[2], ti, similar(x0))
    NonSmoothSDEStep(step1, step2)
end

function (pcl::PreComputeNewtonStep)(vi_sde::SDE_VIO, method::DiscreteTimeStepping{<:ExplicitDriftMethod}, _x0, _x1, _t0, _t1) 
    tracer1 = pcl(vi_sde.sde, method, _x0, _x1, _t0, _t1) # When the dynamics is smooth

    # When there is an impact (k = d = 2)
    
    # @variables x[1:d] y[1:k-1] par[1:length(_par(sde))] t0 t1
    @variables x0 v0 xi vi x1 v1 t0 ti t1 
    @variables par[1:(_par(vi_sde) isa Nothing ? 1 : length(_par(vi_sde)))]
    @syms r(v)
    step_sym_i = eval_driftstep_xI_sym(vi_sde.sde, method, [x0, v0], par, t0, ti)
    step_sym_1 = eval_driftstep_xI_sym(vi_sde.sde, method, [xi, -r(vi)*vi], par, ti, t1)
    _eq = [
        xi - step_sym_i[1],
        vi - step_sym_i[2],
        x1 - step_sym_1[1],
        v1 - step_sym_1[2]
    ]
    _x = [x0, v0, ti, vi]
    J_sym = Symbolics.jacobian(_eq,_x);
    _corr = lu(J_sym)\_eq
    x_new = [_x[i] - _corr[i] for i in 1:4]
    _, x_0i! = build_function(x_new, _x, [x1, v1], par, xi, r, t0, t1, expression = Val{false})
    
    JI_sym = collect(J_sym[1:3, 1:3]) |> lu
    _corr = JI_sym\collect(_eq[1:3])
    x_new = [_x[i] - _corr[i] for i in 1:3]
    _, xI_0i! = build_function(x_new, [x0, ti, vi], [x1, v1], par, xi, r, t0, t1, expression = Val{false})

    detJI_inv = [1/Symbolics.derivative(x1 - substitute(step_sym_1[1], Dict(xi => step_sym_i[1], vi => step_sym_i[2])),x0)]

    detJI⁻¹ = build_function(detJI_inv, [x0, v0], par, xi, r, t0, ti, t1, expression = Val{false})

    
    # _, xII_1! = build_function(step_sym[k:d], x, par, dt, expression = Val{false})
    tracer2 = SymbolicNewtonStepTracer(xI_0i!, x_0i!, detJI⁻¹, similar(_x0,3), similar(_x0,4))
    return (tracer1, tracer2)
end

# function (pcl::PreComputeNewtonStep)(sde::SDE_VIO, method::DiscreteTimeStepping{<:Euler}, x0, x1, _t0, _t1) 
#     step1 = pcl(sde, method, x0, x1, _t0, _t1) # When the dynamics is smooth

#     # When there is an impact (k = d = 2)
    
#     # @variables x[1:d] y[1:k-1] par[1:length(_par(sde))] t0 t1
#     @variables x0 v0 xi vi x1 v1 par[1:length(_par(sde))] t0 ti t1 
#     @syms r(v)
#     step_sym_i = eval_driftstep_xI_sym(sde, method, [x0, v0], par, t0, ti)
#     step_sym_1 = eval_driftstep_xI_sym(sde, method, [xi, -r(vi)*vi], par, ti, t1)
#     _eq = [
#         xi - step_sym_i[1],
#         vi - step_sym_i[2],
#         x1 - step_sym_1[1],
#         v1 - step_sym_1[2]
#     ]
#     _x = [ti, vi, x0, v0]
#     J_sym = Symbolics.jacobian(_eq,_x);
#     _corr = lu(J_sym)\_eq
#     x_new = [_x[i] - _corr[i] for i in 1:4]
#     _, x_0i! = build_function(x_new, _x, [x1, v1], par, xi, r, t0, t1, expression = Val{false})
    
#     JI_sym = collect(J_sym[1:3, 1:3]) |> lu
#     _corr = JI_sym\collect(_eq[1:3])
#     x_new = [_x[i] - _corr[i] for i in 1:3]
#     _, xI_0i! = build_function(x_new, [ti, vi, x0], [x1, v1], par, xi, r, t0, t1, expression = Val{false})

#     detJI_inv = [1/Symbolics.derivative(x1 - substitute(step_sym_1[1], Dict(xi => step_sym_i[1], vi => step_sym_i[2])),x0)]

#     detJI⁻¹ = build_function(detJI_inv, [x0, v0], par, xi, r, t0, ti, t1, expression = Val{false})

    
#     # _, xII_1! = build_function(step_sym[k:d], x, par, dt, expression = Val{false})
#     step2 = SymbolicNewtonStepTracer(xI_0i!, x_0i!, detJI⁻¹, similar(x0,3), similar(x0,4))
#     SymbolicNonSmoothNewtonStepTracer((step1, step2))
# end

##
# driftstep.jl
# TODO: change to Ref(...) for time placeholders in SDEStep
# TODO: get rid of t::Number implementations
# TODO: 
function compute_missing_states_driftstep!(step::SDEStep{d,k,m,sdeT,DiscreteTimeStepping{TDrift,TDiff}}; max_iter = 100, atol = sqrt(eps()), kwargs...) where {d,k,m,sdeT,TDrift, TDiff}
    i = 1
    x_change = 2atol
    while x_change > atol && i < max_iter
        iterate_xI0!(step)
        x_change = norm(step.x0[j] - x for (j,x) in enumerate(step.steptracer.tempI))
        update_drift_x!(step)
        i = i + 1
    end
    # println("iterations: $(i-1)")
    nothing
end
function compute_initial_states_driftstep!(step::SDEStep{d,k,m,sdeT,DiscreteTimeStepping{TDrift,TDiff}}; max_iter = 100, atol = sqrt(eps()), kwargs...) where {d,k,m,sdeT,TDrift, TDiff}
    i = 1
    x_change = 2atol
    while x_change > atol && i < max_iter
        iterate_x0!(step)
        x_change = norm(step.x0[j] - x for (j,x) in enumerate(step.steptracer.temp))
        step.x0 .= step.steptracer.temp
        i = i + 1
    end
    # println("iterations: $(i-1)")
    nothing
end

# _eval_driftstep!
# Near impact: for higher order methods try stepping and if the intermediate step exceeds the impact barrier then recompute it with impact.