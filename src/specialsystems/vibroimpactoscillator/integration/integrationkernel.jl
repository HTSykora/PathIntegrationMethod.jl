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
    if sdestep.Q_aux[] && ID_dyn == 2
        sdestep[ID_dyn].x1[2] = sdestep[ID_dyn].xi2[2]
        # println("ti adjusted")
        sdestep[ID_dyn].ti[] = _t1(sdestep[ID_dyn])
    end
    transitionprobability(sdestep[ID_dyn],x1)
end
function detJ_correction(fx,sdestep::NonSmoothSDEStep{d,k,m})  where {dk,d,k,m}
    fx*get_detJinv(sdestep[sdestep.ID_dyn[]])
end

function detJ_correction(fx,sdestep::NonSmoothSDEStep{d,k,m,sdeT})  where {sdeT<:SDE_VIO,dk,d,k,m}
    return fx * get_detJinv(sdestep[sdestep.ID_dyn[]])
end