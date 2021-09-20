# scalar problem
#TODO
function (IK::IntegrationKernel{sdeT})(vals,x₀) where sdeT<:AbstractSDE{1,1}
    # basefun_vals_safe!(IK.xs.itp,vals,IK.xs,x₀)
    # fx = _tp(IK.sde,IK.pdgrid.xs[1][IK.idx₁...],get_t1(IK), x₀,get_t0(IK), method = IK.method);
    # vals .*= fx
    # vals
end


# multivariate smooth problem
function (IK::IntegrationKernel{sdeT})(vals,x) where sdeT<:AbstractSDE{d,k,m} where {d,k,m}
    update_relevant_states!(IK,x)
    compute_missing_states_driftstep!(IK.sdestep)

    # TODO if m ≠ 1 and d≠k: figure out a rework
    fx = transitionprobability(IK.sdestep,x)

    ## Get interpolation values
    basefun_vals_safe!(IK)
    ## Summarise
    fill_vals!(vals,IK,fx)

    vals
end

# Utility functions:
function update_relevant_states!(IK::IntegrationKernel{sdeT},x::Number) where sdeT<:AbstractSDE{d,d,m} where {d,m}
    IK.sdestep.x0[d] = x
end
function update_relevant_states!(IK::IntegrationKernel{sdeT},x::AbstractVector) where sdeT<:AbstractSDE{d,k,m} where {d,k,m}
    for (i,j) in enumerate(k:d)
        IK.sdestep.x0[j] = x[i]
    end
end

function basefun_vals_safe!(IK::IntegrationKernel{sdeT}) where sdeT<:AbstractSDE{d} where d
    basefun_vals_safe!.(IK.temp.itpVs, IK.pdf.axes,IK.sdestep.x0)
end


get_tempval(str::AbstractVector, i) = str[i]
function fill_vals!(vals::AbstractArray{T,d}, IK::IntegrationKernel{sdeT}, fx) where {T,sdeT<:AbstractSDE{d}} where d
    for idxs in IK.temp.idx_it
        vals[idxs...] = fx * reduce_tempprod(zip(IK.temp.itpVs,idx)...)
        # prod(IK.temp.itpVs[i][idx[i]] for (i,idx) in enumerate(idxs))
    end
    vals
end