struct MeanFirstPassageTime{dynT,TT,cT,IKT,MXT, kwargT}
    sdestep::dynT
    T::TT # Interpolated mean first passage time function T(x)
    condition::cT # condition(x) return true if we are in the state we want to reach
    IK::IKT
    MFPTMX::MXT
    kwargs::kwargT
end

struct MFPTIntegrationKernel{kd,dynT,cT,diT,x0T,x1T,TT,tempT,kwargT}
    sdestep::dynT
    condition::cT
    discreteintegrator::diT
    x0::x0T
    x1::x1T
    T::TT # Interpolated mean first passage time function T(x)
    temp::tempT
    kwargs::kwargT
end
dense_idx_it(IK::MFPTIntegrationKernel) = BI_product(eachindex.(IK.T.axes)...)
_idx_it(IK::MFPTIntegrationKernel) = _idx_it(IK.temp)
_val_it(IK::MFPTIntegrationKernel) = _val_it(IK.temp)

# Fix x0

fill_MFPTMX(MX::Transpose, IK::MFPTIntegrationKernel; kwargs...) = fill_MFPTMX(MX.parent, IK; kwargs...)

function fill_MFPTMX(MX, IK::MFPTIntegrationKernel; smart_integration = true, kwargs...)
    inv_idxs = InvertedIndex(Vector{Int}(undef,0))
    for (i,idx) in enumerate(dense_idx_it(IK))
        update_IK_state_x0_by_idx!(IK, idx)
        
        if IK.condition(IK.x0)
            push!(inv_idxs.skip,i)
            continue  # continue if in region we want to reach
        end

        update_dyn_state_x0!(IK)
        eval_driftstep!(IK)
        update_IK_state_x1_from_sdestep(IK)
        if smart_integration
            rescale_discreteintegrator!(IK; IK.kwargs...)
        end

        get_IK_weights!(IK; IK.kwargs...)
        fill_to_stepMX!(MX,IK,i; IK.kwargs...)
    end
    return inv_idxs, MX[inv_idxs,inv_idxs]
end

function update_IK_state_x0_by_idx!(IK::MFPTIntegrationKernel{kd, dyn}, idx) where {dyn<:AbstractSDEStep{d}} where {kd, d}
    for i in 1:d
        IK.x0[i] = getindex(IK.T.axes[i],idx[i])
    end
end
function update_IK_state_x1_from_sdestep!(IK::MFPTIntegrationKernel{kd, dyn})
    IK.x1 .= IK.sdestep.x1
end
function update_dyn_state_x0!(IK::MFPTIntegrationKernel)
    update_dyn_state_x0!(IK.sdestep,IK.x0)
end

eval_driftstep!(IK::MFPTIntegrationKernel) = eval_driftstep!(IK::sdestep)

function rescale_discreteintegrator!(IK::MFPTIntegrationKernel{1,dyn}; int_limit_thickness_multiplier = 6, kwargs...) where dyn <:SDEStep{d,d,m} where {d,m}
    σ = sqrt(_Δt(IK.sdestep)*get_g(IK.sdestep.sde)(d, IK.sdestep.x0,_par(IK.sdestep),_t0(IK.sdestep))^2)
    mn = min(IK.T.axes[d][end], max(IK.T.axes[d][1],IK.sdestep.x1[d] - int_limit_thickness_multiplier*σ))
    mx = max(IK.T.axes[d][1],min(IK.T.axes[d][end],IK.sdestep.x1[d] + int_limit_thickness_multiplier*σ))
    rescale_to_limits!(IK.discreteintegrator, mn, mx)
end


function get_IK_weights!(IK::MFPTIntegrationKernel; kwargs...)
    IK.discreteintegrator(IK, IK.temp.itpM; kwargs...)
end

function (IK::MFPTIntegrationKernel{dk})(vals,x) where dk
    fx = transitionprobability(sdestep,x)
    # * if m ≠ 1 and d ≠ k: figure out a rework
    set_oldvals_tozero!(vals, IK)
    update_relevant_states!(IK,x)

    if isapprox(_fx,zero(_fx), atol=1e-8) || IK.condition(x)
        all_zero!(vals, IK)
    else
        basefun_vals_safe!(IK)
        fill_vals!(vals,IK,fx,_idx_it(IK), _val_it(IK))
    end
end

# function update_relevant_states!(IK::MFPTIntegrationKernel{dk,sdeT},x) where sdeT<:SDEStep{d,d,m} where {dk,d,m}
#     @inbounds IK.x1[d] = x
# end
function update_relevant_states!(IK::MFPTIntegrationKernel{dk,sdeT},x) where sdeT<:SDEStep{d,k,m} where {dk,d,k,m}
    for (i,j) in enumerate(k:d)
        @inbounds IK.x1[j] = x[i]
    end
end

function basefun_vals_safe!(IK::MFPTIntegrationKernel{dk,sdeT}) where sdeT<:AbstractSDEStep{d} where {dk,d}
    for (it,ax,x1) in zip(IK.temp.itpVs,IK.T.axes,IK.x1)
        basefun_vals_safe!(it,ax,x1; IK.kwargs...)
    end
    nothing
end



