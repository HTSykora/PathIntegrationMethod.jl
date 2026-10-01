function IntegrationKernel(sdestep, f, discreteintegrator, ts, pdf, ikt, kwargs = nothing)
    x1 = similar_to_x1(sdestep)
    kd = getintegration_dimensions(discreteintegrator)

    IntegrationKernel{kd, typeof(sdestep), typeof(x1), typeof(discreteintegrator), typeof(f), typeof(pdf), typeof(ts), typeof(ikt), typeof(kwargs)}(sdestep,x1, f, discreteintegrator, ts, pdf, ikt, kwargs)
end

IK_temp(itpVs, itpM, kernel = GenericRowKernel()) = IK_temp(BI_product(_eachindex.(itpVs)...), BI_product(_val.(itpVs)...), itpVs, itpM, kernel)

# Row kernel used with the discrete integrator `di`
Q_generic_row_kernel(di) = !(di isa AbstractDiscreteIntegratorMethod)
row_kernel(pdf::InterpolatedFunction{T,N,<:SparseInterpolationType}, n_q) where {T,N} = SparseAccumulator(zeros(Bool, length(pdf.p)), Int[])
function row_kernel(pdf::InterpolatedFunction{T,N,<:DenseInterpolationType}, n_q) where {T,N}
    Ns = size(pdf.p)
    Bs = ntuple(j -> zeros(T, Ns[j], n_q), N)
    KR = N ≥ 3 ? zeros(T, prod(Ns[2:end]), n_q) : nothing
    DenseTensorKernel(Bs, zeros(T, n_q), zeros(T, Ns[1], n_q), KR)
end

# Evaluating the integrals

get_IK_weights!(IK::IntegrationKernel; kwargs...) = get_IK_weights!(IK.temp.kernel, IK; kwargs...)
function get_IK_weights!(::GenericRowKernel, IK::IntegrationKernel; kwargs...)
# function get_IK_weights!(IK::IntegrationKernel{1}; integ_limits = first(IK.int_axes), kwargs...)
    IK.discreteintegrator(IK, IK.temp.itpM; kwargs...)
    # quadgk!(IK, IK.temp.itpM, integ_limits...; cleanup_quadgk_keywords(;kwargs...)...)
end

# Sparse interpolations: accumulate w * fx * (interpolation weights) only on the stencil of each quadrature node
function get_IK_weights!(spa::SparseAccumulator, IK::IntegrationKernel; kwargs...)
    itpM = IK.temp.itpM
    reset_row!(spa, itpM)
    di = IK.discreteintegrator
    if di.Q_integrate[]
        LI = LinearIndices(itpM)
        for (w,x) in zip(di.w,di.x)
            fx = kernel_value!(IK, x)
            if !iszero(fx)
                basefun_vals_safe!(IK)
                accumulate_vals!(spa, itpM, LI, fx, w, _idx_it(IK), _val_it(IK))
            end
        end
    end
    nothing
end
function reset_row!(spa::SparseAccumulator, itpM)
    for j in spa.touched
        itpM[j] = zero(eltype(itpM))
        spa.mark[j] = false
    end
    empty!(spa.touched)
    nothing
end
function accumulate_vals!(spa::SparseAccumulator, itpM, LI, fx, w, idx_it, val_it)
    for (idx, val) in zip(idx_it, val_it)
        v = fx * prod(val) * w # the same operations as in `fill_vals!` and in the DiscreteIntegrator
        if !iszero(v)
            j = LI[idx...]
            if !spa.mark[j]
                spa.mark[j] = true
                push!(spa.touched, j)
            end
            itpM[j] += v
        end
    end
    nothing
end

# Dense interpolations: Φ = Σ_q c_q ⊗_j Bs[j][:,q], evaluated as a matrix product
function get_IK_weights!(dk::DenseTensorKernel, IK::IntegrationKernel; kwargs...)
    itpM = IK.temp.itpM
    di = IK.discreteintegrator
    if !di.Q_integrate[]
        fill!(itpM, zero(eltype(itpM)))
        return nothing
    end
    for (q, (w,x)) in enumerate(zip(di.w,di.x))
        fx = kernel_value!(IK, x)
        if iszero(fx)
            dk.c[q] = zero(eltype(dk.c))
            for B in dk.Bs
                fill!(view(B, :, q), zero(eltype(B)))
            end
        else
            dk.c[q] = w * fx
            x0 = get_sdestep_x0(IK)
            map_axes((B, ax, i) -> basefun_vals_safe!(view(B, :, q), ax, x0[i]; IK.kwargs...), dk.Bs, IK.pdf.axes)
        end
    end
    contract_row!(itpM, dk)
    nothing
end
contract_row!(itpM::AbstractVector, dk::DenseTensorKernel) = mul!(itpM, dk.Bs[1], dk.c)
function contract_row!(itpM::AbstractMatrix, dk::DenseTensorKernel)
    dk.B1c .= dk.Bs[1] .* transpose(dk.c)
    mul!(itpM, dk.B1c, transpose(dk.Bs[2]))
end
function contract_row!(itpM::AbstractArray, dk::DenseTensorKernel)
    dk.B1c .= dk.Bs[1] .* transpose(dk.c)
    khatrirao!(dk.KR, Base.tail(dk.Bs))
    mul!(reshape(itpM, size(itpM, 1), :), dk.B1c, transpose(dk.KR))
end
# KR[:,q] = kron(Bs[end][:,q], …, Bs[1][:,q]), i.e. the index of Bs[1] runs fastest (column-major order)
function khatrirao!(KR, Bs)
    CI = CartesianIndices(map(B -> axes(B, 1), Bs))
    for q in axes(KR, 2)
        for (l, I) in enumerate(CI)
            KR[l, q] = prod(ntuple(j -> Bs[j][I[j], q], length(Bs)))
        end
    end
    KR
end

# multivariate smooth problem
function (IK::IntegrationKernel{dk,sdeT})(vals,x) where sdeT<:AbstractSDEStep{d,k,m} where {dk,d,k,m}
    # * if m ≠ 1 and d ≠ k: figure out a rework
    set_oldvals_tozero!(vals, IK)
    fx = kernel_value!(IK, x)
    if iszero(fx)
        all_zero!(vals, IK)
    else
        basefun_vals_safe!(IK)
        fill_vals!(vals,IK,fx,_idx_it(IK), _val_it(IK))
    end

    vals
end

# Transitional PDF with the Jacobian correction at the integration variable `x` (zero if negligible)
function kernel_value!(IK::IntegrationKernel, x)
    _fx = _getTPDF(IK.sdestep, x, IK.x1)
    if isapprox(_fx,zero(_fx), atol=1e-8)
        return zero(_fx)
    end
    detJ_correction(_fx,IK.sdestep)
end

function _getTPDF(sdestep, x, x1)
    update_relevant_states!(sdestep,x)
    compute_missing_states_driftstep!(sdestep)
    transitionprobability(sdestep,x1)
end

# Utility functions:
function detJ_correction(fx,sdestep::SDEStep{d,1,m})where {d,m}
    fx
end
function detJ_correction(fx,sdestep::SDEStep{d,k,m}) where {k,d,m}
    fx*get_detJinv(sdestep)
end

function update_relevant_states!(IK::IntegrationKernel{dk,sdeT},x) where sdeT<:SDEStep{d,k,m} where {dk,d,k,m,N}
        update_relevant_states!(IK.sdestep, x)
end

function update_relevant_states!(sdestep::sdeT,x::Number) where sdeT<:SDEStep{d,d,m} where {dk,d,m}
    @inbounds sdestep.x0[d] = x
end
# the noisy coordinates k:d (several integration variables: an SVector)
function update_relevant_states!(sdestep::SDEStep{d,k,m}, x::SVector) where {d,k,m}
    for (i,j) in enumerate(k:d)
        @inbounds sdestep.x0[j] = x[i]
    end
end
function update_relevant_states!(sdestep::sdeT,x::Vararg{Any,N}) where sdeT<:SDEStep{d,k,m} where {dk,d,k,m,N}
    for (i,j) in enumerate(k:d)
        @inbounds sdestep.x0[j] = x[i]
    end
end

function set_oldvals_tozero!(vals, IK::IntegrationKernel)
end
function set_oldvals_tozero!(vals, IK::IntegrationKernel{kd,sdeT,x1T,xT,fT,pdfT}) where {kd,sdeT,x1T,xT,fT,pdfT<:InterpolatedFunction{T,N, itp_type}} where {T,N,itp_type <: SparseInterpolationType}
    @inbounds for idx in _idx_it(IK)
        vals[idx...] = zero(eltype(vals))
    end
end
function all_zero!(vals::AbstractArray{<:T}, IK::IntegrationKernel) where T
    vals .= zero(T)
    nothing
end
function all_zero!(vals::AbstractArray{<:T}, IK::IntegrationKernel{kd,sdeT,x1T,xT,fT,pdfT}) where {kd,sdeT,x1T,xT,fT,pdfT<:InterpolatedFunction{T,N, itp_type}} where {T,N,itp_type <: SparseInterpolationType}
end

function basefun_vals_safe!(IK::IntegrationKernel{dk,sdeT}) where sdeT<:AbstractSDEStep{d} where {dk,d}
    x0 = get_sdestep_x0(IK)
    map_axes((it, ax, i) -> basefun_vals_safe!(it, ax, x0[i]; IK.kwargs...), IK.temp.itpVs, IK.pdf.axes)
    nothing
end
# f(a[i], axes[i], i) for each axis: the loop over the tuples is unrolled, as the axes (and their buffers) can have different
# types (e.g. a QuinticAxis and a ChebyshevAxis); iterating over them would be type unstable
map_axes(f, a::NTuple{N,Any}, axes::NTuple{N,Any}) where N = (map(f, a, axes, ntuple(identity, Val(N))); nothing)

function get_sdestep_x0(IK::IntegrationKernel)
    IK.sdestep.x0
end

_idx_it(IK::IntegrationKernel) = _idx_it(IK.temp)
_idx_it(IKT::IK_temp) = IKT.idx_it
dense_idx_it(IK::IntegrationKernel) = BI_product(eachindex.(IK.pdf.axes)...)

_val_it(IK::IntegrationKernel) = _val_it(IK.temp)
_val_it(IKT::IK_temp) = IKT.val_it

# get_tempval(str::AbstractVector, i) = str[i]
function fill_vals!(vals::AbstractArray{T,d}, IK::IntegrationKernel{dk,sdeT}, fx, idx_it, val_it;) where {T,sdeT<:AbstractSDEStep{d}} where {dk,d}
    for (idx, val) in zip(idx_it, val_it)
        vals[idx...] = fx * prod(val)# reduce_tempprod(zip(IK.temp.itpVs,idx)...)
        # prod(IK.temp.itpVs[i][idx[i]] for (i,idx) in enumerate(idxs))
    end
    nothing
    #vals
end

