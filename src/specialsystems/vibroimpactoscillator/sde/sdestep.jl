similar_to_x1(sdestep::NonSmoothSDEStep, args...) = similar(first(sdestep.sdesteps).x1, args...)

function set_t0t1!(sdestep::NonSmoothSDEStep,t0,t1) 
    for _sdestep in sdestep.sdesteps
        set_t0t1!(_sdestep,t0,t1)
    end
end

function NonSmoothSDEStep(sde::sdeT, sdesteps::Vararg{SDEStep, N}; kwargs...) where {sdeT, N}
    @assert reduce(&, Q_compatible(sdesteps[1], sdesteps[i]) for i in 2:N) "Incompatible SDE steps provided"
    d,k,m = get_dkm(sdesteps[1])
    ID_dyn = Ref(1)
    ID_aux = nothing
    Q_aux = Ref(false)
    NonSmoothSDEStep{d,k,m,sdeT, N, typeof(sdesteps),  typeof(ID_dyn), typeof(ID_aux), typeof(Q_aux)}(sde, sdesteps, ID_dyn, ID_aux, Q_aux)
end
Base.getindex(sdestep::NonSmoothSDEStep,idx...) = sdestep.sdesteps[idx...]
Base.size(::NonSmoothSDEStep{d,k,m,sdeT,n}) where {d,k,m,sdeT,n} = n
function set_wall_ID!(sde::SDE_VIO,ID)
    sde.ID[] = ID
end
function set_wall_ID!(sdestep::SDEStep{d,k,m,sdeT},ID) where {d,k,m,sdeT<:SDE_VIO}
    set_wall_ID!(sdestep.sde,ID)
end
function set_wall_ID!(sdestep::NonSmoothSDEStep{d,k,m,sdeT},ID) where {d,k,m,sdeT<:SDE_VIO}
    set_wall_ID!(sdestep[2],ID)
end
get_wall_ID(sde::SDE_VIO) = sde.ID[]
get_wall_ID(sdestep::SDEStep{d,k,m,sdeT}) where {d,k,m,sdeT<:SDE_VIO} = get_wall_ID(sdestep.sde)
get_wall_ID(sdestep::NonSmoothSDEStep{d,k,m,sdeT}) where {d,k,m,sdeT<:SDE_VIO} = get_wall_ID(sdestep[2])
get_wall(sde::SDE_VIO) = sde.wall[sde.ID[]]
get_wall(sdestep::SDEStep{d,k,m,sdeT}) where {d,k,m,sdeT<:SDE_VIO} = get_wall(sdestep.sde)
get_wall(sdestep::NonSmoothSDEStep{d,k,m,sdeT}) where {d,k,m,sdeT<:SDE_VIO} = get_wall(sdestep[2])

function SDEStep(sde::sdeT, method::methodT, x0,x1, t0, t1; precomputelevel::pclT = PreComputeNewtonStep(), kwargs...) where {sdeT<:SDE_VIO, methodT <: DiscreteTimeSteppingMethod, pclT <: PreComputeLevel} where {d,k,m}
    
    _method = DiscreteTimeStepping(sde, method)
    steptracers = precomputelevel(sde,_method, x0,x1, t0, t1)
    ti = init_ti(t0,t1)
        
    step1 = SDEStep{2,2,1,typeof(sde.sde),typeof(_method),typeof(steptracers[1]),typeof(x0),typeof(x1),typeof(t0), Nothing, Nothing, typeof(x0)}(sde.sde, _method, similar(x0), similar(x1), t0, t1, steptracers[1], nothing, nothing, similar(x0))
    step2 = SDEStep{2,2,1,sdeT,typeof(_method),typeof(steptracers[2]),typeof(x0),typeof(x1),typeof(t0), typeof(ti), typeof(x0), typeof(x0)}(sde, _method, similar(x0), similar(x1), t0, t1, steptracers[2], ti, similar(x0), similar(x0))
    preset_xi_vals!(step2)
    NonSmoothSDEStep(sde, step1, step2; kwargs...)
end
function init_ti(t0,t1)
    _t0 = get_val(t0)
    Ref(_t0 + (get_val(t1) - _t0)/2)
    # Ref(zero(get_val(t0)))
end
get_val(t::Number) = t
get_val(t::Base.RefValue) = t[]

preset_xi_vals!(step::SDEStep) = step
function preset_xi_vals!(step2::SDEStep{2,2,1,sdeT}) where sdeT<:SDE_VIO{1}
    step2.xi[1] = step2.sde.wall[1].pos
    step2.xi2[1] = step2.xi[1]
    step2
end

function substitute_w_to_rdiff(exprs, W, r, v)
    [expand_derivatives(substitute(expr,Dict(r(v)=>W(v)))) for expr in exprs] |> collect
end
function build_inplace_diffsubstituted_function(expr, W, r, v, args...; kwargs...)
    new_expr = substitute_w_to_rdiff(expr, W, r, v)
    _, f! = build_function(new_expr, args...; expression = Val{false})
    return f!
end
function build_inplace_function(expr, args...; kwargs...)
    _, f! = build_function(expr, args...; expression = Val{false})
    return f!
end

function build_diffsubstituted_function(expr, W, r, v, args...; kwargs...)
    new_expr = expand_derivatives(substitute(expr,Dict(r(v)=>W(v))))
    build_function(new_expr, args...; expression = Val{false})
end
function (pcl::PreComputeNewtonStep)(vi_sde::SDE_VIO, method::DiscreteTimeStepping{<:ExplicitDriftMethod}, _x0, _x1, _t0, _t1) 
    tracer1 = pcl(vi_sde.sde, method, _x0, _x1, _t0, _t1) # When the dynamics is smooth

    # When there is an impact (k = d = 2)
    @variables x0 v0 xi vi x1 v1 t0 ti t1 
    @variables par[1:(_par(vi_sde) isa Nothing ? 1 : length(_par(vi_sde)))]
    @syms r(v)
    step_sym_i = eval_driftstep_xI_sym(vi_sde.sde, method, [x0, v0], par, t0, ti)
    step_sym_1s = Tuple(eval_driftstep_xI_sym(vi_sde.sde, method, [xi, -W(vi)*vi], par, ti, t1) for W in vi_sde.wall)
    _eqs = Tuple([
        xi - step_sym_i[1],
        vi - step_sym_i[2],
        x1 - step_sym_1[1],
        v1 - step_sym_1[2]
    ] for step_sym_1 in step_sym_1s)
    _x = (x0, v0, ti, vi)
    J_syms =Tuple([Symbolics.derivative(eq,x) for eq in _eq, x in _x] for _eq in _eqs)#Symbolics.jacobian(_eq,_x);
    _corrs = Tuple(lu(J_sym)\_eq for( J_sym,_eq) in zip(J_syms,_eqs))
    x_news = Tuple([_x[i] - _corr[i] for i in 1:4] for _corr in _corrs)
    x_0i! = Tuple(build_inplace_function(x_new, [x0, v0], [x1, v1], [xi, vi], par, t0, t1, ti) for x_new in x_news)

    _x = (x0,  ti, vi)
    _eq13s = Tuple(collect(_eq[1:3]) for _eq in _eqs)
    JI_syms = Tuple([Symbolics.derivative(eq,x) for eq in _eq13, x in _x] for _eq13 in _eq13s)
    
    _corrs = Tuple(JI_sym\_eq13 for (JI_sym,_eq13) in zip(JI_syms,_eq13s))
    x_news = Tuple([_x[i] - _corr[i] for i in 1:3] for _corr in _corrs)
    xI_0i! = Tuple(build_inplace_function(x_new, [x0, v0], [x1, v1], [xi, vi], par, t0, t1, ti) for x_new in x_news)


    dtidx0 = -Symbolics.derivative(step_sym_i[1], x0)/Symbolics.derivative(step_sym_i[1], ti)
    fullstep_sym_1s = Tuple(substitute(step_sym_1[1], Dict(xi => step_sym_i[1], vi => step_sym_i[2])) for step_sym_1 in step_sym_1s)
    detJI_invs = Tuple(-1/(Symbolics.derivative(step_sym_1, x0) + Symbolics.derivative(step_sym_1, ti)*dtidx0) for step_sym_1 in fullstep_sym_1s)
    
    detJI⁻¹ = Tuple(build_function(detJI_inv, [x0, v0], [x1, v1], [xi, vi], par, t0, t1, ti; expression = Val{false}) for detJI_inv in detJI_invs)

    _xi = [v0, vi, v1]
    # v0 is the velocity just before the impact
    step_sym_ai = eval_driftstep_xI_sym(vi_sde.sde, method, [xi, -r(v0)*v0], par, t0, t1)
    _eq_i = [
        x1 - step_sym_ai[1],
        v1 - step_sym_ai[2],
        vi + r(v0)*v0
    ]
    J_i_sym = Symbolics.jacobian(_eq_i,_xi)
    _corr_i = lu(J_i_sym) \ _eq_i
    x_new_i = [_xi[i] - _corr_i[i] for i in 1:3]
    x_i! = Tuple(build_inplace_diffsubstituted_function(x_new_i, W, r, v0, [v0, vi, v1], x1, xi, par, t0, t1) for W in vi_sde.wall)
    # _, xII_1! = build_function(step_sym[k:d], x, par, dt, expression = Val{false})
    tracer2 = VIO_SymbolicNewtonImpactStepTracer(xI_0i!, x_0i!, detJI⁻¹, similar(_x0,3), similar(_x0,4), x_i!, similar(_x0,3), similar(_x0,3))
    return (tracer1, tracer2)
end

function apply_correction_to_xI0!(step::SDEStep{d,k,m, sdeT, methodT,tracerT}; kwargs...) where {d, k,m, sdeT, methodT,tracerT<:VIO_SymbolicNewtonImpactStepTracer}
    step.steptracer.xI_0![get_wall_ID(step)](step.steptracer.tempI, step.x0, step.x1, step.xi, _par(step), _t0(step),_t1(step), _ti(step))
end
function apply_correction_to_x0!(step::SDEStep{d,k,m, sdeT, methodT,tracerT}; kwargs...) where {d, k,m, sdeT, methodT,tracerT<:VIO_SymbolicNewtonImpactStepTracer}
    step.steptracer.x_0![get_wall_ID(step)](step.steptracer.temp, step.x0, step.x1, step.xi, _par(step), _t0(step),_t1(step), _ti(step))
end

# To update the 
function _compute_velocities_to_impact!(temp, v_f!, vs, x1, xi, par, t0, t1; max_iter = 100, atol = 1.5e-8, kwargs...)
    i = 1
    x_change = 2atol
    while x_change > atol && i < max_iter
        v_f!(temp, vs, x1, xi, par, t0, t1; kwargs...)
        x_change = norm(vs[j] - temp[j] for j in eachindex(vs))
        for j in eachindex(vs)
            vs[j] = temp[j]
        end
        i = i + 1
    end
    nothing
end
function compute_velocities_to_impact!(step::SDEStep{2,2,1,<:SDE_VIO,DiscreteTimeStepping{TDrift,TDiff}}; kwargs...) where {TDrift, TDiff}
    # Compute max v0, v1 just before impact happens
    _compute_velocities_to_impact!(step.steptracer.vitemp,step.steptracer.v_toimpact![get_wall_ID(step)], step.steptracer.v_i, step.x1[1], step.xi[1], _par(step), _t0(step), _t1(step); kwargs...)
    abs(step.steptracer.v_i[1])<1.5e-8, step.steptracer.v_i
    # nothing
end

function get_detJinv(step::SDEStep{d,k,m, sdeT, methodT,tracerT}) where {d, k,m, sdeT<:SDE_VIO, methodT,tracerT<:VIO_SymbolicNewtonImpactStepTracer}
    step.steptracer.detJI_inv[get_wall_ID(step)](step.x0, step.x1, step.xi, _par(step), _t0(step),_t1(step), _ti(step)) |> abs
end
