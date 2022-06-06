
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
function NonSmoothSDEStep(sde::sdeT, sdesteps::Vararg{SDEStep, N}; Q_switch = Ref(false), ID = Ref(0), kwargs...) where {sdeT, N}
    @assert reduce(&, Q_compatible(sdesteps[1], sdesteps[i]) for i in 2:N) "Incompatible SDE steps provided"
    d,k,m = get_dkm(sdesteps[1])
    NonSmoothSDEStep{d,k,m,sdeT, typeof(sdesteps), typeof(Q_switch), typeof(ID)}(sde, sdesteps, Q_switch, ID)
end

function SDEStep(sde::sdeT, method::methodT, x0,x1, t0, t1; precomputelevel::pclT = PreComputeNewtonStep(), kwargs...) where {sdeT<:SDE_VIO, methodT <: DiscreteTimeSteppingMethod, pclT <: PreComputeLevel} where {d,k,m}
    
    _method = DiscreteTimeStepping(sde, method)
    steptracers = precomputelevel(sde,_method, x0,x1, t0, t1)
    if t0 isa Base.RefValue
        ti = Ref(zero(typeof(t0[])));
    else
        ti = Ref(zero(typeof(t0)));
    end
        
    step1 = SDEStep{2,2,1,sdeT,typeof(_method),typeof(steptracers[1]),typeof(x0),typeof(x1),typeof(t0), Nothing, Nothing, Nothing}(sde, _method, similar(x0), similar(x1), t0, t1, steptracers[1], nothing, nothing, nothing)
    step2 = SDEStep{2,2,1,sdeT,typeof(_method),typeof(steptracers[2]),typeof(x0),typeof(x1),typeof(t0), typeof(ti), typeof(x0), typeof(x0)}(sde, _method, similar(x0), similar(x1), t0, t1, steptracers[2], ti, similar(x0), similar(x0))
    NonSmoothSDEStep(sde, step1, step2; kwargs...)
end

function substitute_w_to_rdiff(exprs, W, r, v)
    [expand_derivatives(substitute(expr,Dict(r(v)=>W(v)))) for expr in exprs] |> collect
end
function build_inplace_diffsubstituted_function(expr, W, r, v, args...; kwargs...)
    new_expr = substitute_w_to_rdiff(expr, W, r, v)
    _, f! = build_function(new_expr, args...; expression = Val{false})
    return f!
end
function (pcl::PreComputeNewtonStep)(vi_sde::SDE_VIO, method::DiscreteTimeStepping{<:ExplicitDriftMethod}, _x0, _x1, _t0, _t1) 
    tracer1 = pcl(vi_sde.sde, method, _x0, _x1, _t0, _t1) # When the dynamics is smooth

    # When there is an impact (k = d = 2)
    # TODO: substitute wall functions instead of r(vi)
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
    _x = (x0, v0, ti, vi)
    J_sym = [Symbolics.derivative(eq,x) for eq in _eq, x in _x]#Symbolics.jacobian(_eq,_x);
    _corr = lu(J_sym)\_eq
    x_new = [_x[i] - _corr[i] for i in 1:4]
    x_0i! = Tuple(build_inplace_diffsubstituted_function(x_new, W, r, vi, [x0, v0], [x1, v1], [xi, vi], par, t0, t1, ti) for W in vi_sde.wall)

    _x = (x0,  ti, vi)
    _eq13 = collect(_eq[1:3])
    JI_sym = [Symbolics.derivative(eq,x) for eq in _eq13, x in _x]
    # JI_sym = collect(J_sym[1:3, 1:3]) |> lu
    _corr = JI_sym\_eq13
    x_new = [_x[i] - _corr[i] for i in 1:3]
    xI_0i! = Tuple(build_inplace_diffsubstituted_function(x_new, W, r, vi, [x0, v0], [x1, v1], [xi, vi], par, t0, t1, ti) for W in vi_sde.wall)

    detJI_inv = [1/Symbolics.derivative(x1 - substitute(step_sym_1[1], Dict(xi => step_sym_i[1], vi => step_sym_i[2])),x0)]

    detJI⁻¹ =Tuple(build_inplace_diffsubstituted_function(detJI_inv, W, r, vi, [x0, v0], [x1, v1], [xi, vi], par, t0, t1, ti) for W in vi_sde.wall)
    # detJI⁻¹ = build_function(detJI_inv, [x0, v0], [xi, vi], par, r, t0, t1, ti, expression = Val{false})

    _xi = [v0, v1]
    step_sym_bi = eval_driftstep_xI_sym(vi_sde.sde, method, [xi, v0], par, t0, t1)
    step_sym_ai = eval_driftstep_xI_sym(vi_sde.sde, method, [xi, -r(v0)*v0], par, t0, t1)
    _eq_bi = [
        x1 - step_sym_bi[1],
        v1 - step_sym_bi[2]
    ]
    _eq_ai = [
        x1 - step_sym_ai[1],
        v1 - step_sym_ai[2]
    ]
    J_bi_sym = Symbolics.jacobian(_eq_bi,_xi)
    J_ai_sym = Symbolics.jacobian(_eq_ai,_xi)
    _corr_bi = lu(J_bi_sym) \ _eq_bi
    _corr_ai = lu(J_ai_sym) \ _eq_ai
    x_new_bi = [_xi[i] - _corr_bi[i] for i in 1:2]
    x_new_ai = [_xi[i] - _corr_ai[i] for i in 1:2]
    #TODO: figure out the necessary arguments
    xI_0bi! = Tuple(build_inplace_diffsubstituted_function(x_new_bi, W, r, vi, [x0, v0], [x1, v1], [xi, vi], par, t0, t1, ti) for W in vi_sde.wall)
    xI_0ai! = Tuple(build_inplace_diffsubstituted_function(x_new_ai, W, r, vi, [x0, v0], [x1, v1], [xi, vi], par, t0, t1, ti) for W in vi_sde.wall)
    # _, xII_1! = build_function(step_sym[k:d], x, par, dt, expression = Val{false})
    tracer2 = VIO_SymbolicNewtonImpactStepTracer(xI_0i!, x_0i!, detJI⁻¹, similar(_x0,3), similar(_x0,4), xI_0bi!, xI_0ai!, similar(_x0,2))
    return (tracer1, tracer2)
end


function iterate_xI0!(step::SDEStep{d,k,m, sdeT, methodT,tracerT}; wallID = 1) where {d, k,m, sdeT, methodT,tracerT<:VIO_SymbolicNewtonImpactStepTracer}
    step.steptracer.xI_0![wallID](step.steptracer.tempI, step.x0, step.x1, step.xi, _par(step), _t0(step),_t1(step), _ti(step))
end

# function update_relevant_states!(IK::IntegrationKernel{dk,sdeT},v0::Number) where sdeT<:SDE_VIO where {dk}
#     IK.sdestep.sdesteps[1].x0[2] = v0
#     update_impact_states_if_necessary!(IK.sdestep, v0, IK.sdestep.sde.wall...; IK.kwargs...)
# end
# function update_impact_states_if_necessary!(nonsmoothsdestep::NonSmoothSDEStep, v0, wall; impact_range_modifier = get_impact_range_modifier(nonsmoothsdestep.sdesteps[2].method), kwargs...)
#     _x0 = nonsmoothsdestep.sdesteps[2].x1[1] + impact_range_modifier*wall(v0)*v0*_Δt(sdestep) 
#     if _x0 < wall.pos
#         nonsmoothsdestep.Q_switch[] = true
#         nonsmoothsdestep.ID[] = 1
#         nonsmoothsdestep.sdesteps[2].x0[2] = v0
#     end
# end
# function update_impact_states_if_necessary!(nonsmoothsdestep::NonSmoothSDEStep, v0, wall1, wall2; impact_range_modifier = get_impact_range_modifier(nonsmoothsdestep.sdesteps[2].method), kwargs...)
#     _x0 = nonsmoothsdestep.sdesteps[2].x1[1] + impact_range_modifier*wall1(v0)*v0*_Δt(sdestep) 
#     if _x0 < wall.pos
#         nonsmoothsdestep.Q_switch[] = true
#         nonsmoothsdestep.ID[] = 1
#         nonsmoothsdestep.sdesteps[2].x0[2] = v0
#         return
#     end
#     _x0 = nonsmoothsdestep.sdesteps[2].x1[1] + impact_range_modifier*wall2(v0)*v0*_Δt(sdestep) 
#     if _x0 > wall.pos
#         nonsmoothsdestep.Q_switch[] = true
#         nonsmoothsdestep.ID[] = 1
#         nonsmoothsdestep.sdesteps[2].x0[2] = v0
#         return
#     end
# end
# get_impact_range_modifier(::DiscreteTimeStepping{<:Euler}) = 1.
# get_impact_range_modifier(::DiscreteTimeStepping{<:RungeKutta}) = 1.1


##

# driftstep.jl
# TODO: change to Ref(...) for time placeholders in SDEStep
# TODO: get rid of t::Number implementations
# TODO: 
function compute_missing_states_driftstep!(step::NonSmoothSDEStep{d,k,m,sdeT}; kwargs...) where {d,k,m,sdeT<:SDE_VIO}
    if step.Q_switch[]
        compute_missing_states_driftstep!(step.sdesteps[2],update_impact_vio_x!, wallID = step.ID[])
    else
        compute_missing_states_driftstep!(step.sdesteps[1])
    end
end
# function compute_missing_states_vio_driftstep!(step::SDEStep{d,k,m,sdeT,DiscreteTimeStepping{TDrift,TDiff}}; max_iter = 100, atol = sqrt(eps()), kwargs...) where {d,k,m,sdeT,TDrift, TDiff}
#     i = 1
#     x_change = 2atol
#     while x_change > atol && i < max_iter
#         iterate_xI0!(step)
#         x_change = norm(step.x0[j] - step.steptracer.tempI[j] for j in 1:(k-1))
        
#         update_impact_vio_x!(step, kwargs...)
#         i = i + 1
#     end
#     # println("iterations: $(i-1)")
#     nothing
# end

function update_impact_vio_x!(step::SDEStep{2,2,1, sdeT, methodT,tracerT,x0T,x1T,tT, tiT, xiT, xiT}; wallID = 1) where {sdeT<:SDE_VIO,methodT<:DiscreteTimeStepping{TDrift},tracerT,x0T,x1T,tT, tiT, xiT} where {TDrift}
    # for i in 1:k-1
    #     step.x0[i] = step.steptracer.tempI[i]
    # end
    step.x0[1] = step.steptracer.tempI[1]
    step.ti[] = step.steptracer.tempI[2]
    step.xi[2] = step.steptracer.tempI[3] # vi
    step.xi2[2] = - step.sde.wall[wallID](step.xi[2])*step.xi[2]
    update_x1_kd!(step)
end

function update_x1_kd!(step::SDEStep{2,2,1, sdeT, methodT,tracerT,x0T,x1T,tT}) where {sdeT<:SDE_VIO,methodT<:DiscreteTimeStepping{TDrift},tracerT,x0T,x1T,tT} where {TDrift<:Euler}
    i = 2; # i in k:d
    step.x1[i]  = step.xi2[i] + get_f(step.sde)(i,step.xi2,_par(step),_ti(step)) * _Δti1(step)
end
function update_x1_kd!(step::SDEStep{d,k,m, sdeT, methodT,tracerT,x0T,x1T,tT}) where {d,k,m,sdeT<:SDE_VIO,methodT<:DiscreteTimeStepping{TDrift},tracerT,x0T,x1T,tT} where {TDrift<:RungeKutta{ord}} where ord
    _eval_driftstep!(step,step.xi2,_Δti1(step))
    fill_to_x1!(step,step.xi2,k)
end
