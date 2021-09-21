function IntegrationKernel(sdestep, f, int_axes, ts, pdf, ikt, kwargs = nothing)
    IntegrationKernel{length(int_axes), typeof(sdestep), typeof(int_axes), typeof(f), typeof(pdf), typeof(ts), typeof(ikt), typeof(kwargs)}(sdestep, f, int_axes, ts, pdf, ikt, kwargs)
end

# Evaluating the integrals
@inline function cleanup_quadgk_keywords(;σ_init = nothing, μ_init = nothing, kwargs...)
    (;kwargs...)
end
function get_IK_weights!(IK::IntegrationKernel{1}; integ_limits = (first(IK.int_axes)[1],first(IK.int_axes)[end]), kwargs...)
    quadgk!(IK, IK.temp.itpM, integ_limits...; cleanup_quadgk_keywords(;kwargs...)...)
end

# multivariate smooth problem
function (IK::IntegrationKernel{dk,sdeT})(vals,x) where sdeT<:SDEStep{d,k,m} where {dk,d,k,m}
    update_relevant_states!(IK,x)
    compute_missing_states_driftstep!(IK.sdestep)

    # * if m ≠ 1 and d ≠ k: figure out a rework
    fx = transitionprobability(IK.sdestep,x)

    ## Get interpolation values
    basefun_vals_safe!(IK)

    ## Summarise
    fill_vals!(vals,IK,fx,_idx_it(IK))

    vals
end

# Utility functions:
function update_relevant_states!(IK::IntegrationKernel{dk,sdeT},x::Number) where sdeT<:SDEStep{d,d,m} where {dk,d,m}
    IK.sdestep.x0[d] = x
end
function update_relevant_states!(IK::IntegrationKernel{dk,sdeT},x::AbstractVector) where sdeT<:SDEStep{d,k,m} where {dk,d,k,m}
    for (i,j) in enumerate(k:d)
        IK.sdestep.x0[j] = x[i]
    end
end

function basefun_vals_safe!(IK::IntegrationKernel{dk,sdeT}) where sdeT<:SDEStep{d} where {dk,d}
    
    basefun_vals_safe!.(IK.temp.itpVs, IK.pdf.axes,IK.sdestep.x0)
    nothing
end

_idx_it(IK::IntegrationKernel) = _idx_it(IK.temp)
_idx_it(IKT::IK_temp) = IKT.idx_it

get_tempval(str::AbstractVector, i) = str[i]
function fill_vals!(vals::AbstractArray{T,d}, IK::IntegrationKernel{dk,sdeT}, fx, idx_it;) where {T,sdeT<:SDEStep{d}} where {dk,d}
    for idx in idx_it
        vals[idx...] = fx * reduce_tempprod(zip(IK.temp.itpVs,idx)...)
        # prod(IK.temp.itpVs[i][idx[i]] for (i,idx) in enumerate(idxs))
    end
    nothing
    #vals
end

