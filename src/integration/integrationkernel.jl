function IntegrationKernel(sdestep, f, discreteintegrator, ts, pdf, ikt, kwargs = nothing)
    x1 = similar(sdestep.x1)
    if discreteintegrator isa QuadGKIntegrator
        kd = 1
    else
        kd = length(first(discreteintegrator.x))
    end
    IntegrationKernel{kd, typeof(sdestep), typeof(x1), typeof(discreteintegrator), typeof(f), typeof(pdf), typeof(ts), typeof(ikt), typeof(kwargs)}(sdestep,x1, f, discreteintegrator, ts, pdf, ikt, kwargs)
end

# Evaluating the integrals

function get_IK_weights!(IK::IntegrationKernel{1})
# function get_IK_weights!(IK::IntegrationKernel{1}; integ_limits = first(IK.int_axes), kwargs...)
    IK.discreteintegrator(IK, IK.temp.itpM)
    # quadgk!(IK, IK.temp.itpM, integ_limits...; cleanup_quadgk_keywords(;kwargs...)...)
end

# multivariate smooth problem
function (IK::IntegrationKernel{dk,sdeT})(vals,x) where sdeT<:SDEStep{d,k,m} where {dk,d,k,m}
    update_relevant_states!(IK,x)
    compute_missing_states_driftstep!(IK.sdestep)

    # * if m ≠ 1 and d ≠ k: figure out a rework
    fx = transitionprobability(IK.sdestep,IK.x1)

    ## Get interpolation values
    set_oldvals_tozero!(vals, IK)
    if isapprox(fx,zero(fx),atol=1e-8)
        all_zero!(vals, IK)
    else
        basefun_vals_safe!(IK)
        fill_vals!(vals,IK,fx,_idx_it(IK), _val_it(IK))
    end

    vals
end

# Utility functions:
function update_relevant_states!(IK::IntegrationKernel{dk,sdeT},x::Number) where sdeT<:SDEStep{d,d,m} where {dk,d,m}
    IK.sdestep.x0[d] = x
end
function update_relevant_states!(IK::IntegrationKernel{dk,sdeT},x::Vararg{Any,N}) where sdeT<:SDEStep{d,k,m} where {dk,d,k,m,N}
    for (i,j) in enumerate(k:d)
        IK.sdestep.x0[j] = x[i]
    end
end

function set_oldvals_tozero!(vals, IK::IntegrationKernel)
end
function set_oldvals_tozero!(vals, IK::IntegrationKernel{kd,sdeT,x1T,xT,fT,pdfT}) where {kd,sdeT,x1T,xT,fT,pdfT<:InterpolatedFunction{T,N, itp_type}} where {T,N,itp_type <: SparseInterpolationType}
    for idx in _idx_it(IK)
        vals[idx...] = zero(eltype(vals))
    end
end
function all_zero!(vals::AbstractArray{<:T}, IK::IntegrationKernel) where T
    vals .= zero(T)
    nothing
end
function all_zero!(vals::AbstractArray{<:T}, IK::IntegrationKernel{kd,sdeT,x1T,xT,fT,pdfT}) where {kd,sdeT,x1T,xT,fT,pdfT<:InterpolatedFunction{T,N, itp_type}} where {T,N,itp_type <: SparseInterpolationType}
end

function basefun_vals_safe!(IK::IntegrationKernel{dk,sdeT}) where sdeT<:SDEStep{d} where {dk,d}
    for (it,ax,x0) in zip(IK.temp.itpVs,IK.pdf.axes,IK.sdestep.x0)
        basefun_vals_safe!(it,ax,x0)
    end
    nothing
end

_idx_it(IK::IntegrationKernel) = _idx_it(IK.temp)
_idx_it(IKT::IK_temp) = IKT.idx_it
dense_idx_it(IK::IntegrationKernel) = BI_product(eachindex.(IK.pdf.axes)...)

_val_it(IK::IntegrationKernel) = _val_it(IK.temp)
_val_it(IKT::IK_temp) = IKT.val_it

# get_tempval(str::AbstractVector, i) = str[i]
function fill_vals!(vals::AbstractArray{T,d}, IK::IntegrationKernel{dk,sdeT}, fx, idx_it, val_it;) where {T,sdeT<:SDEStep{d}} where {dk,d}
    for (idx, val) in zip(idx_it, val_it)
        vals[idx...] = fx * prod(val)# reduce_tempprod(zip(IK.temp.itpVs,idx)...)
        # prod(IK.temp.itpVs[i][idx[i]] for (i,idx) in enumerate(idxs))
    end
    nothing
    #vals
end

