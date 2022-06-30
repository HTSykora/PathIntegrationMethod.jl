
###################
(w::Wall{<:Function})(v) = w.r(abs(v))
(w::Wall{<:Number})(v) = w.r
Wall(r) = Wall(r,0.);

# sde/sde.jl
get_dkm(sde::SDE_VIO) = (2,2,1) # (d,k,m)

function osc_f1(u,p,t)
    u[2]
end
is_v_compatible(v, ID) = !xor(signbit(v), ID == 1)
is_x_inimpactzone(x, xi, ID; wall_tol = 1.5e-8, kwargs...) = 
    isapprox(x, xi, atol = wall_tol) ? false : xor(ID == 1, x<xi)

function SDE_VIO(f::fT,g::gT, wall::Union{Tuple{wT1},Tuple{wT1,wT2}}, par=nothing) where {fT<:Function,gT<:Function, wT1<:Wall, wT2<:Wall}
    sde = SDE((osc_f1,f), g, par);
    ID = Ref(1)
    SDE_VIO{length(wall),typeof(sde), typeof(wall), typeof(ID)}(sde, wall, ID)
end
SDE_VIO(f::fT,g::gT, w::Wall, par = nothing; kwargs...) where {fT<:Function,gT<:Function} = SDE_VIO(f, g, (w,), par; kwargs...)
get_Q_cv(::Tuple{wT1}) where wT1<:Wall = (signbit,)
get_Q_cv(::Tuple{wT1,wT2}) where {wT1<:Wall, wT2<:Wall} = (signbit, !signbit)

# sde/sdestep.jl
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
    ID_aux = get_ID_aux(sde)
    Q_aux = Ref(false)
    NonSmoothSDEStep{d,k,m,sdeT, N, typeof(sdesteps),  typeof(ID_dyn), typeof(ID_aux), typeof(Q_aux)}(sde, sdesteps, ID_dyn, ID_aux, Q_aux)
end
Base.getindex(sdestep::NonSmoothSDEStep,idx...) = sdestep.sdesteps[idx...]
Base.size(::NonSmoothSDEStep{d,k,m,sdeT,n}) where {d,k,m,sdeT,n} = n
get_ID_aux(sde::AbstractSDE) = Ref(1)
get_ID_aux(sde::SDE_VIO) = sde.ID



function SDEStep(sde::sdeT, method::methodT, x0,x1, t0, t1; precomputelevel::pclT = PreComputeNewtonStep(), kwargs...) where {sdeT<:SDE_VIO, methodT <: DiscreteTimeSteppingMethod, pclT <: PreComputeLevel} where {d,k,m}
    
    _method = DiscreteTimeStepping(sde, method)
    steptracers = precomputelevel(sde,_method, x0,x1, t0, t1)
    if t0 isa Base.RefValue
        ti = Ref(zero(typeof(t0[])));
    else
        ti = Ref(zero(typeof(t0)));
    end
        
    step1 = SDEStep{2,2,1,typeof(sde.sde),typeof(_method),typeof(steptracers[1]),typeof(x0),typeof(x1),typeof(t0), Nothing, Nothing, typeof(x0)}(sde.sde, _method, similar(x0), similar(x1), t0, t1, steptracers[1], nothing, nothing, similar(x0))
    step2 = SDEStep{2,2,1,sdeT,typeof(_method),typeof(steptracers[2]),typeof(x0),typeof(x1),typeof(t0), typeof(ti), typeof(x0), typeof(x0)}(sde, _method, similar(x0), similar(x1), t0, t1, steptracers[2], ti, similar(x0), similar(x0))
    preset_xi_vals!(step2)
    NonSmoothSDEStep(sde, step1, step2; kwargs...)
end
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
    # JI_sym = collect(J_sym[1:3, 1:3]) |> lu
    _corrs = Tuple(JI_sym\_eq13 for (JI_sym,_eq13) in zip(JI_syms,_eq13s))
    x_news = Tuple([_x[i] - _corr[i] for i in 1:3] for _corr in _corrs)
    xI_0i! = Tuple(build_inplace_function(x_new, [x0, v0], [x1, v1], [xi, vi], par, t0, t1, ti) for x_new in x_news)

    # TODO:
    # step_sym_1_Ws = Tuple(Tuple(substitute(expr,Dict(r(vi)=>W(vi),xi => step_sym_i[1], vi => step_sym_i[2])) for expr in step_sym_1) for W in vi_sde.wall)
    detJI_invs = Tuple(1/Symbolics.derivative(x1 - substitute(step_sym_1[1], Dict(xi => step_sym_i[1], vi => step_sym_i[2])), x0) for step_sym_1 in step_sym_1s)


    detJI⁻¹ = Tuple(build_function(detJI_inv, [x0, v0], [x1, v1], [xi, vi], par, t0, t1, ti; expression = Val{false}) for detJI_inv in detJI_invs)
    # detJI⁻¹ = build_function(detJI_inv, [x0, v0], [xi, vi], par, r, t0, t1, ti, expression = Val{false})

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

function apply_correction_to_xI0!(step::SDEStep{d,k,m, sdeT, methodT,tracerT}; wallID = 1, kwargs...) where {d, k,m, sdeT, methodT,tracerT<:VIO_SymbolicNewtonImpactStepTracer}
    step.steptracer.xI_0![wallID](step.steptracer.tempI, step.x0, step.x1, step.xi, _par(step), _t0(step),_t1(step), _ti(step))
end
function apply_correction_to_x0!(step::SDEStep{d,k,m, sdeT, methodT,tracerT}; wallID = 1, kwargs...) where {d, k,m, sdeT, methodT,tracerT<:VIO_SymbolicNewtonImpactStepTracer}
    step.steptracer.x_0![wallID](step.steptracer.temp, step.x0, step.x1, step.xi, _par(step), _t0(step),_t1(step), _ti(step))
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
function compute_velocities_to_impact!(step::SDEStep{2,2,1,<:SDE_VIO,DiscreteTimeStepping{TDrift,TDiff}}, wallID; kwargs...) where {TDrift, TDiff}
    # Compute max v0, v1 just before impact happens
    _compute_velocities_to_impact!(step.steptracer.vitemp,step.steptracer.v_toimpact![wallID], step.steptracer.v_i, step.x1[1], step.xi[1], _par(step), _t0(step), _t1(step); kwargs...)
    nothing
end

function get_detJinv(step::SDEStep{d,k,m, sdeT, methodT,tracerT}, wallID) where {d, k,m, sdeT<:SDE_VIO, methodT,tracerT<:VIO_SymbolicNewtonImpactStepTracer}
    step.steptracer.detJI_inv[wallID](step.x0, step.x1, step.xi, _par(step), _t0(step),_t1(step), _ti(step)) |> abs
end

##

# driftstep.jl
compute_missing_states_driftstep!(step::SDEStep{d,k,m,sdeT}; kwargs...) where {d,k,m,sdeT<:SDE_VIO} = compute_missing_states_driftstep!(step, update_impact_vio_xI!; kwargs...)
# compute_missing_states_driftstep!(step.sdesteps[2],update_impact_vio_xI!, wallID = step.ID_aux[]; kwargs...)
# function compute_missing_states_driftstep!(step::NonSmoothSDEStep{d,k,m,sdeT}; kwargs...) where {d,k,m,sdeT<:SDE_VIO}
#     if step.ID_dyn[] == 1
#         compute_missing_states_driftstep!(step.sdesteps[1]; kwargs...)
#     else 
#         compute_missing_states_driftstep!(step.sdesteps[2],update_impact_vio_xI!, wallID = step.ID_aux[]; kwargs...)
#     end
# end
# function update_impact_vio_xi_atwall!(step::SDEStep{2,2,1, sdeT, methodT,tracerT,x0T,x1T,tT, tiT, xiT, xiT}, wallID; kwargs...) where {sdeT<:SDE_VIO,methodT<:DiscreteTimeStepping{TDrift},tracerT,x0T,x1T,tT, tiT, xiT} where {TDrift}
#     xi = step.sde.wall[wallID].pos
#     step.xi[1] = xi
#     step.xi2[1] = xi
# end

function update_impact_vio_xI!(step::SDEStep{2,2,1, sdeT, methodT,tracerT,x0T,x1T,tT, tiT, xiT, xiT}; wallID = 1) where {sdeT<:SDE_VIO,methodT<:DiscreteTimeStepping{TDrift},tracerT,x0T,x1T,tT, tiT, xiT} where {TDrift}

    step.x0[1] = step.steptracer.tempI[1]

    step.ti[] = step.steptracer.tempI[2]

    step.xi[2] = step.steptracer.tempI[3] # vi
    step.xi2[2] = - step.sde.wall[wallID](step.xi[2])*step.xi[2]

    update_x1_kd!(step)
end
function update_impact_vio_x!(step::SDEStep{2,2,1, sdeT, methodT,tracerT,x0T,x1T,tT, tiT, xiT, xiT}; wallID = 1, kwargs...) where {sdeT<:SDE_VIO,methodT<:DiscreteTimeStepping{TDrift},tracerT,x0T,x1T,tT, tiT, xiT} where {TDrift}

    step.x0[1] = step.steptracer.temp[1]
    step.x0[2] = step.steptracer.temp[2]

    step.ti[] = step.steptracer.temp[3]

    step.xi[2] = step.steptracer.temp[4] # vi
    step.xi2[2] = - step.sde.wall[wallID](step.xi[2])*step.xi[2]
end


function compute_initial_states_driftstep!(step2::SDEStep{d,k,m,sdeT}; wallID = 1, kwargs...) where {d,k,m,sdeT<:SDE_VIO}
    compute_initial_states_driftstep!(step2, update_impact_vio_x!; wallID = wallID, kwargs...)
end

function update_x1_kd!(step::SDEStep{2,2,1, sdeT, methodT,tracerT,x0T,x1T,tT}) where {sdeT<:SDE_VIO,methodT<:DiscreteTimeStepping{TDrift},tracerT,x0T,x1T,tT} where {TDrift<:Euler}
    i = 2; # i in k:d
    step.x1[i]  = step.xi2[i] + get_f(step.sde)(i,step.xi2,_par(step),_ti(step)) * _Δti1(step)
end
function update_x1_kd!(step::SDEStep{d,k,m, sdeT, methodT,tracerT,x0T,x1T,tT}) where {d,k,m,sdeT<:SDE_VIO,methodT<:DiscreteTimeStepping{TDrift},tracerT,x0T,x1T,tT} where {TDrift<:RungeKutta{ord}} where ord
    _eval_driftstep!(step,step.xi2,_Δti1(step))
    fill_to_x1!(step,step.xi2,k)
end

##
# compute_stepMX.jl
function _v_beforeimpact_f(v1, v0, r)
    v1 + v0*r(v0)
end
function _v_new(v1, v2, f1, f2)
    v1 - f1 * (v2 - v1) / (f2-f1)
end
function get_v_beforeimpact(v_afterimpact, wall::Wall{<:Number}; kwargs...)
    - v_afterimpact / wall.r
end
function get_v_beforeimpact(v_afterimpact, wall::Wall; v_maxiter = 100, v0_neighbourhood_half = 0.025, v_atol = 1.5e-8, kwargs...)
    v_new = - v_afterimpact / wall(v_afterimpact)
    v_min = v_new*(one(v_new) - v0_neighbourhood_half)
    v_max = v_new*(one(v_new) + v0_neighbourhood_half)

    f_at_v_min = _v_beforeimpact_f(v_afterimpact, v_min, wall)
    f_at_v_max = _v_beforeimpact_f(v_afterimpact, v_max, wall)
    
    xchange = 2v_atol
    v = v_new
    i = 1
    while xchange > v_atol && i < v_maxiter
        v_new = _v_new(v_min, v_max, f_at_v_min, f_at_v_max)
        xchange = abs(v - v_new)
        v = v_new
        f_at_vnew =  _v_beforeimpact_f(v_afterimpact, v_new, wall)
        if sign(f_at_vnew) == sign(f_at_v_min)
            v_min = v_new
            f_at_v_min = _v_beforeimpact_f(v_afterimpact, v_min, wall)
        else
            v_max = v_new
            f_at_v_max = _v_beforeimpact_f(v_afterimpact, v_max, wall)
        end
        i = i + 1
    end
    v
end

function get_and_set_potential_wallID!(IK::IntegrationKernel{1,dyn}) where dyn <:NonSmoothSDEStep{2,2,1, sdeT} where {sdeT<:SDE_VIO}
    get_and_set_potential_wallID!(IK.sdestep, IK.x1[1])
end
function get_and_set_potential_wallID!(sdestep::NonSmoothSDEStep{2,2,1,sdeT}, x1) where sdeT <: SDE_VIO{1}
    # sdestep.ID_aux[] = 1
end
function get_and_set_potential_wallID!(sdestep::NonSmoothSDEStep{2,2,1,sdeT}, x1) where sdeT <: SDE_VIO{2}
    walls = sdestep.sde.wall
    ID = abs(walls[1].pos - x1) < abs(walls[2].pos - x1) ? 1 : 2
    sdestep.ID_aux[] = ID
    W_pos = walls[ID].pos
    # update_impact_vio_xi_atwall!(sdestep[2], ID)

    sdestep[2].xi[1] = W_pos
    sdestep[2].xi2[1] = W_pos
end
function modify_x1_if_on_wall!(IK::IntegrationKernel{1,dyn}; wall_tol = 1.5e-8, kwargs... ) where dyn <:NonSmoothSDEStep{2,2,1, sdeT} where {sdeT<:SDE_VIO}
    # sqrt(eps()) ≈ 1.49e-8
    ID = IK.sdestep.ID_aux[];
    wall = IK.sdestep.sde.wall[ID[]];
    x1,v1 = IK.x1
    if (isapprox(x1, wall.pos, atol = wall_tol) && !(is_v_compatible(v1, ID)))
        IK.sdestep.sdesteps[1].xi2[1] = x1
        IK.sdestep.sdesteps[1].xi2[2] = v1
        IK.x1[2] = get_v_beforeimpact(v1, wall; kwargs...)
        IK.sdestep.Q_aux[] = true
    end
end

function update_IK_state_x1!(IK::IntegrationKernel{1,dyn}, idx) where dyn <:NonSmoothSDEStep{2,2,1, sdeT} where {sdeT<:SDE_VIO}
    d=2
    IK.sdestep.Q_aux[] = false

    for i in 1:d
        IK.x1[i] = getindex(IK.pdf.axes[i],idx[i])
    end

    get_and_set_potential_wallID!(IK)
    modify_x1_if_on_wall!(IK)
end
function update_dyn_state_x1!(IK::IntegrationKernel{1,dyn}, idx) where dyn <:NonSmoothSDEStep{d,k,m, sdeT} where {d,k,m,sdeT}
    for sdestep in IK.sdestep.sdesteps
        update_dyn_state_x1!(sdestep,IK.x1)
    end
end


function Q_check_impact(step, walls, ID)
    if ID == 1
        return step.x0[1] < walls[1].pos
    else
        return step.x0[1] > walls[2].pos
    end
end

function rescale_discreteintegrator!(IK::IntegrationKernel{1,dyn}; int_limit_thickness_multiplier = 6, kwargs...) where dyn <:NonSmoothSDEStep{2,2,1, sdeT} where {sdeT<:SDE_VIO}
    
    IK.discreteintegrator[1].Q_integrate[] = false
    IK.discreteintegrator[2].Q_integrate[] = false
    
    sdestep = IK.sdestep
    step1 = sdestep.sdesteps[1]
    step2 = sdestep.sdesteps[2]
    compute_velocities_to_impact!(step2, sdestep.ID_aux[])
    s_σ = sqrt(_Δt(step1)*get_g(step1.sde)(2, step1.x1,_par(step1),_t0(step1))^2) * int_limit_thickness_multiplier


    if abs(step2.steptracer.v_i[3] - IK.x1[2]) ≤ s_σ && !(sdestep.Q_aux[]) # v1 after impact vs. v1 investigated
        # mixed: impact and non-impact
        wallID = sdestep.ID_aux[]
        if wallID == 1
            # without impact
            mx = min(step2.steptracer.v_i[2],IK.pdf.axes[2][end])
            step1.x1[2] = step1.x1[2] - s_σ
            compute_initial_states_driftstep!(step1; IK.kwargs...)
            mn = max(step1.x0[2],IK.pdf.axes[2][1])
            rescale_to_limits!(IK.discreteintegrator[1], mn, mx)
            step1.x1[2] = step2.x1[2]
            
            # with impact
            mx = min(step2.steptracer.v_i[1],IK.pdf.axes[2][end])
            step2.x1[2] = step2.x1[2] + s_σ
            compute_initial_states_driftstep!(step2; wallID = wallID, IK.kwargs...)
            mn = max(step2.x0[2],IK.pdf.axes[2][1])
            rescale_to_limits!(IK.discreteintegrator[2], mn, mx)
            step2.x1[2] = step1.x1[2]
        elseif wallID == 2
            # without impact
            mn = max(step2.steptracer.v_i[2],IK.pdf.axes[2][1])
            step1.x1[2] = step1.x1[2] + s_σ
            compute_initial_states_driftstep!(step1; IK.kwargs...)
            mx = min(step1.x0[2],IK.pdf.axes[2][end])
            rescale_to_limits!(IK.discreteintegrator[1], mn, mx)
            step1.x1[2] = step2.x1[2]
            
            # with impact
            mn = max(step2.steptracer.v_i[1],IK.pdf.axes[2][1])
            step2.x1[2] = step2.x1[2] - s_σ
            compute_initial_states_driftstep!(step2; wallID = wallID, IK.kwargs...)
            mx = min(step2.x0[2],IK.pdf.axes[2][end])
            rescale_to_limits!(IK.discreteintegrator[2], mn, mx)
            step2.x1[2] = step1.x1[2]
        end
    else # check if there is a clean impact region or not?
        compute_initial_states_driftstep!(step1; IK.kwargs...)
        if Q_check_impact(step1, step2.sde.wall, sdestep.ID_aux[])
            compute_initial_states_driftstep!(step2; IK.kwargs...)
            rescale_discreteintegrator!(IK.discreteintegrator[2], step2, IK.pdf; int_limit_thickness_multiplier = int_limit_thickness_multiplier, kwargs...)
        else
            # checkbit = sdestep.Q_aux[] && abs(step1.x1[2])<1.4e-5
            rescale_discreteintegrator!(IK.discreteintegrator[1], step1, IK.pdf; int_limit_thickness_multiplier = int_limit_thickness_multiplier, kwargs...)
            # if checkbit
            #     @show IK.discreteintegrator[1].x[1]
            #     @show IK.discreteintegrator[1].x[end]
            #     @show IK.discreteintegrator[1].Q_integrate
            #     @show IK.discreteintegrator[2].Q_integrate
            # end
        end
    end
end
## 

# function get_IK_weights!(IK::IntegrationKernel{1})
#     # function get_IK_weights!(IK::IntegrationKernel{1}; integ_limits = first(IK.int_axes), kwargs...)
#     IK.discreteintegrator.discreteintegrators(IK, IK.temp.itpM)
#     # quadgk!(IK, IK.temp.itpM, integ_limits...; cleanup_quadgk_keywords(;kwargs...)...)
# end
    
##
# integrationkernel.jl
function get_sdestep_x0(IK::IntegrationKernel{dk,stepT}) where {dk,stepT<:NonSmoothSDEStep}
    IK.sdestep[IK.sdestep.ID_dyn[]].x0
end
function get_IK_weights!(IK::IntegrationKernel{1,stepT}; kwargs...)  where stepT<:NonSmoothSDEStep{d,k,m,sdeT,n} where {d,k,m,sdeT,n}
    # function get_IK_weights!(IK::IntegrationKernel{1}; integ_limits = first(IK.int_axes), kwargs...)
    IK.sdestep.ID_dyn[] = 1
    IK.discreteintegrator[1](IK, IK.temp.itpM; Q_reinit_res = true,  kwargs...)
    for i in 2:n
        IK.sdestep.ID_dyn[] = i
        IK.discreteintegrator[i](IK, IK.temp.itpM;  Q_reinit_res = false)
    end
        # quadgk!(IK, IK.temp.itpM, integ_limits...; cleanup_quadgk_keywords(;kwargs...)...)
end
function _getTPDF(sdestep::NonSmoothSDEStep{d,k,m}, x, x1)  where {dk,d,k,m}
    ID_dyn = sdestep.ID_dyn[]
    update_relevant_states!(sdestep[ID_dyn],x)
    compute_missing_states_driftstep!(sdestep[ID_dyn])
    transitionprobability(sdestep[ID_dyn],x1)
end
function detJ_correction(fx,sdestep::NonSmoothSDEStep{d,k,m})  where {dk,d,k,m}
    fx*get_detJinv(sdestep[sdestep.ID_dyn[]])
end

function detJ_correction(fx,sdestep::NonSmoothSDEStep{d,k,m,sdeT})  where {sdeT<:SDE_VIO,dk,d,k,m}
    if sdestep.Q_aux[]
        return fx * get_detJinv(sdestep[2], sdestep.ID_aux[])
    else
        return fx * get_detJinv(sdestep[sdestep.ID_dyn[]], sdestep.ID_aux[])
    end
end

# discreteintegrator.jl
Base.getindex(di::NonSmoothDiscreteIntegrator, idx...) = di.discreteintegrators[idx...]
Base.size(::NonSmoothDiscreteIntegrator{n,NoDyn}) where {n,NoDyn}= NoDyn
function DiscreteIntegrator(discreteintegrator, sdestep::NonSmoothSDEStep, res_prototype, N::Union{NTuple{1,<:Integer},<:Integer,AbstractArray{<:Integer}}, axes::GA; kwargs...) where GA
    discreteintegrators = Tuple(DiscreteIntegrator(discreteintegrator,res_prototype, N, axes; kwargs...) for _ in sdestep.sdesteps)
    
    NonSmoothDiscreteIntegrator{1,number_of_sdesteps(sdestep),typeof(discreteintegrators)}(discreteintegrators)
end

# transitionprobabilities.jl
# function transitionprobability(step::SDEStep{d,d,m,sdeT,method},x) where {d,m,sdeT<:SDE_VIO, method<:DiscreteTimeStepping{TDrift, TDiff}} where {TDrift, TDiff<:Maruyama}
#     g = get_g(step.sde)
#     r = step.sde.wall[step.sde.ID[]]
#     σ2 = _Δt0i(step) * (g(d, step.x0,_par(step),_t0(step))^ 2)
#     σ2 = σ2 * (r(step.xi[2])^2)
#     σ2 = σ2 + _Δti1(step) * (g(d, step.xi2,_par(step),_ti(step))^ 2) 
#     # detJ_correction = _detJ(step, x)
#     normal1D_σ2(step.x1[d], σ2, x[d]);
# end