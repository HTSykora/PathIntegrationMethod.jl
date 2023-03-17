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
get_and_set_potential_wallID!(::NonSmoothSDEStep{2,2,1,sdeT}, x1) where sdeT <: SDE_VIO{1} = nothing
function get_and_set_potential_wallID!(sdestep::NonSmoothSDEStep{2,2,1,sdeT}, x1) where sdeT <: SDE_VIO{2}
    walls = sdestep.sde.wall
    ID = abs(walls[1].pos - x1) < abs(walls[2].pos - x1) ? 1 : 2
    set_wall_ID!(sdestep,ID)
    W_pos = get_wall(sdestep).pos
    # update_impact_vio_xi_atwall!(sdestep[2], ID)

    sdestep[2].xi[1] = W_pos
    sdestep[2].xi2[1] = W_pos

    # sdestep[2].Q_aux[] = abs(W_pos - x1) < 1e-10
    nothing
end

function update_IK_state_x1_by_idx!(IK::IntegrationKernel{1,dyn}, idx) where dyn <:NonSmoothSDEStep{2,2,1, sdeT} where {sdeT<:SDE_VIO}
    d=2

    for i in 1:d
        IK.x1[i] = getindex(IK.pdf.axes[i],idx[i])
    end

    get_and_set_potential_wallID!(IK)
end
function update_dyn_state_xs!(IK::IntegrationKernel{1,dyn}) where dyn <:NonSmoothSDEStep{d,k,m, sdeT} where {d,k,m,sdeT}
    for sdestep in IK.sdestep.sdesteps
        update_dyn_state_xs!(sdestep,IK.x1)
    end
end
function update_dyn_state_xs!(sdestep::SDEStep{d,k,m,sdeT}, x1) where {d,k,m,sdeT<:SDE_VIO}
    _corr = get_wall_ID(sdestep) == 1 ? 1.5e-10 : -1.5e-10
    sdestep.x1[1] = x1[1] + _corr
    sdestep.x1[2] = x1[2] + _corr
    sdestep.xi2[2] = x1[2] + _corr
    sdestep.xi[2] = -sdestep.xi2[2]/get_wall(sdestep)(sdestep.xi2[2])

    sdestep.x0[1] = x1[1] + _corr
    sdestep.x0[2] = sdestep.xi[2]

    sdestep.ti[] = _t1(sdestep)

end

function Q_check_impact(step, walls, ID)
    if ID == 1
        return step.x0[1] < walls[1].pos
    else
        return step.x0[1] > walls[2].pos
    end
end

function rescale_discreteintegrator!(IK::IntegrationKernel{1,dyn}; impact_int_limit_thickness_multiplier = 6, kwargs...) where dyn <:NonSmoothSDEStep{2,2,1, sdeT} where {sdeT<:SDE_VIO}
    
    IK.discreteintegrator[1].Q_integrate[] = false
    IK.discreteintegrator[2].Q_integrate[] = false
    
    sdestep = IK.sdestep
    step1 = sdestep.sdesteps[1]
    step2 = sdestep.sdesteps[2]

    Q_at_wall, v_i = compute_velocities_to_impact!(step2)
    # v_i == [v_beforeimpact, v_afterimpact, v1]
    IK.sdestep.Q_aux[] = Q_at_wall
    
    if Q_at_wall
        v_i[3] = zero(v_i[3])
        v_i[1] = get_wall_ID(step2) == 1 ?  -1.5e-8 : 1.5e-8
        v_i[2] =  -v_i[1]*get_wall(step2)(v_i[1])
        step2.ti[] = _t1(step2)
    end
    s_σ = sqrt(_Δt(step1)*get_g(step1.sde)(2, step1.x1,_par(step1),_t0(step1))^2) * impact_int_limit_thickness_multiplier
    # TODO: handle impact when v_0 == 0 at wall!

    if abs(v_i[3] - IK.x1[2]) ≤ s_σ # v1 after impact vs. v1 investigated
        # mixed: impact and non-impact
        wallID = get_wall_ID(sdestep)
        if wallID == 1
            # without impact
            mn, mx = get_rescale_limits(step1, IK.pdf; kwargs..., IK.kwargs...)
            mx = min(v_i[2], mx)
            rescale_to_limits!(IK.discreteintegrator[1], mn, mx)
            
            # with impact
            mx = min(v_i[1], IK.pdf.axes[2][end])
            step2.x1[2] = step2.x1[2] + s_σ
            compute_initial_states_driftstep!(step2; IK.kwargs...)
            mn = max(step2.x0[2],IK.pdf.axes[2][1])
            rescale_to_limits!(IK.discreteintegrator[2], mn, mx)
            step2.x1[2] = step1.x1[2]
        elseif wallID == 2
            # without impact
            mn, mx = get_rescale_limits(step1, IK.pdf; kwargs..., IK.kwargs...)
            mn = max(v_i[2], mn)
            rescale_to_limits!(IK.discreteintegrator[1], mn, mx)
            
            # with impact
            mn = max(v_i[1], IK.pdf.axes[2][1])
            step2.x1[2] = step2.x1[2] - s_σ
            compute_initial_states_driftstep!(step2; IK.kwargs...)
            mx = min(step2.x0[2],IK.pdf.axes[2][end])
            rescale_to_limits!(IK.discreteintegrator[2], mn, mx)
            step2.x1[2] = step1.x1[2]
        end
    else # check if there is a clean impact region or not?
        compute_initial_states_driftstep!(step1; IK.kwargs...)
        if Q_check_impact(step1, step2.sde.wall, get_wall_ID(sdestep))
            compute_initial_states_driftstep!(step2; IK.kwargs...)
            rescale_discreteintegrator!(IK.discreteintegrator[2], step2, IK.pdf; kwargs...)
        else
            rescale_discreteintegrator!(IK.discreteintegrator[1], step1, IK.pdf; kwargs...)
            # if sdestep.Q_aux[]
            #     (mn, mx) = get_limits(IK.discreteintegrator[1])
            #     if get_wall_ID(sdestep) == 1
            #         mx = min(zero(mx),mx)
            #     else
            #         mn = max(zero(mn),mn)
            #     end
            #     rescale_to_limits!(IK.discreteintegrator[2], mn, mx)
            # end
        end
    end
    if Q_at_wall
        update_dyn_state_xs!(step2,step1.x1)
    end
end