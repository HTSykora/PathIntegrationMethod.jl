struct MeanFirstPassageTime{dynT,TT,cT,IKT,MXT, kwargT}
    sdestep::dynT
    T::TT # Interpolated mean first passage time function T(x)
    condition::cT # condition(x) return true if we are in the state we want to reach
    IK::IKT
    MFPTMX::MXT
    kwargs::kwargT
end
function MeanFirstPassageTime(sdestep::AbstractSDEStep{d,k,m}, condition::cT, axes::Vararg{Any,d};
    di_N = 31, discreteintegrator = defaultdiscreteintegrator(sdestep.sde,di_N),
    mfptMXtype = nothing, sparse_tol = 1e-6, kwargs...) where {cT<:Function, d,k,m}

    if mfptMXtype isa StepMatrixRepresentation
        _mfptMXtype = mfptMXtype
    else
        _mfptMXtype = get_stepMXtype(sdestep.sde, get_val_itp_type(axes); sparse_tol = sparse_tol, kwargs...)
    end

    T = InterpolatedFunction(axes...; kwargs...)

    itpVs = Tuple(zero(axis.temp) for axis in axes);
    ikt = IK_temp(BI_product(_eachindex.(itpVs)...), # = idx_it
                        BI_product(_val.(itpVs)...), # = val_it
                        itpVs, # = itpVs
                        zero(T.p)) # = itpM
    di = DiscreteIntegrator(discreteintegrator, sdestep, T.p, axes[k:end]...; kwargs...)
    IK = MFPTIntegrationKernel(sdestep, condition, di, T, ikt, kwargs)
    return compute_MFPT(IK; mfptMXtype = _mfptMXtype)
end

function compute_MFPT(mfptIK; mfptMXtype = DenseMX(), kwargs...)
    _mfptMX = initialize_stepMX(eltype(mfptIK.T.p), length(mfptIK.T),mfptMXtype)

    inv_idxs, mfptMX, mfptV = fill_MFPTMX!(_mfptMX, mfptIK; kwargs...)
    _Tmx = (I - mfptMX') \ mfptV
    # _Tmx = (I - mfptMX') \ fill(_Δt(mfptIK.sdestep), size(mfptMX,1))
    # return inv_idxs, mfptMX
    mfptIK.T.p[inv_idxs] .= vec(_Tmx)
    mfptIK
    # stepMX
end

abstract type AbstractBarrierCondition <: Function end
struct SingleBarrierCondition{i,cT,xT} <: AbstractBarrierCondition
    condition::cT
    x0::xT
end

function SingleBarrierCondition(i::Int, x0; condition = >)
    _cond = condition(x0)
    SingleBarrierCondition{i,typeof(_cond), typeof(x0)}(_cond, x0)
end

function (sc::SingleBarrierCondition{i})(x) where i
    sc.condition(x[i])
end

struct DoubleBarrierCondition{i,c1T,c2T,xT} <: AbstractBarrierCondition
    condition1::c1T
    condition2::c2T
    xs::xT
end

function DoubleBarrierCondition(i::Int, x0, x1)
    _cond1 = <(x0)
    _cond2 = >(x1)
    xs = (x0,x1)
    DoubleBarrierCondition{i,typeof(_cond1),typeof(_cond2), typeof(xs)}(_cond1, _cond2, xs)
end
function DoubleBarrierCondition(i::Int, x0)
    DoubleBarrierCondition(i, -x0, x0)
end

function (dc::DoubleBarrierCondition{i})(x) where i
    dc.condition1(x[i]) || dc.condition2(x[i])
end


abstract type AbstractIntegrationKernel end
struct MFPTIntegrationKernel{kd,dynT,cT,diT,x0T,x1T,TT,tempT,kwargT} <: AbstractIntegrationKernel
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

function MFPTIntegrationKernel(sdestep, condition, discreteintegrator, Tfunc, ik_temp, kwargs = nothing)
    x0 = similar_to_x1(sdestep)
    x1 = similar_to_x1(sdestep)
    kd = getintegration_dimensions(discreteintegrator)
    
    MFPTIntegrationKernel{kd, typeof(sdestep), typeof(condition), typeof(discreteintegrator), typeof(x0), typeof(x1), typeof(Tfunc), typeof(ik_temp), typeof(kwargs)}(sdestep, condition, discreteintegrator, x0, x1, Tfunc, ik_temp, kwargs)
end

function get_sdestep_x0(IK::AbstractIntegrationKernel)
    IK.sdestep.x0
end
function get_sdestep_x1(IK::AbstractIntegrationKernel)
    IK.sdestep.x1
end


fill_MFPTMX!(MX::Transpose, IK::MFPTIntegrationKernel; kwargs...) = fill_MFPTMX!(MX.parent, IK; kwargs...)

function fill_MFPTMX!(MX, IK::MFPTIntegrationKernel; smart_integration = true, kwargs...)
    inv_idxs = InvertedIndex(Vector{Int}(undef,0))
    mfptV = Vector{typeof(_Δt(IK.sdestep))}(undef, 0); sizehint!(mfptV, size(MX,1))
    for (i,idx) in enumerate(dense_idx_it(IK))
        update_IK_state_x0_by_idx!(IK, idx)
        if IK.condition(IK.x0)
            push!(inv_idxs.skip,i)
            # println("Skipped: $(i)")
            continue  # continue if in region we want to reach
        end
        # println("Not skipped: $(i)")
        update_dyn_state_x0!(IK)
        eval_driftstep!(IK)
        update_IK_state_x1_from_sdestep!(IK)
        push_Δt!(mfptV, IK)
        if smart_integration
            rescale_discreteintegrator!(IK; restrict_limit_to_interpolationgrid = false, IK.kwargs...)
        end

        get_IK_weights!(IK; IK.kwargs...)
        fill_to_stepMX!(MX,IK,i; IK.kwargs...)
    end
    return inv_idxs, MX[inv_idxs,inv_idxs], mfptV
end

function update_IK_state_x0_by_idx!(IK::MFPTIntegrationKernel{kd, dyn}, idx) where {dyn<:AbstractSDEStep{d}} where {kd, d}
    for i in 1:d
        IK.x0[i] = getindex(IK.T.axes[i],idx[i])
    end
end
function update_IK_state_x1_from_sdestep!(IK::MFPTIntegrationKernel{kd, dyn}) where {kd,dyn}
    IK.x1 .= IK.sdestep.x1
end
function update_dyn_state_x0!(IK::MFPTIntegrationKernel)
    update_dyn_state_x0!(IK.sdestep,IK.x0)
end

eval_driftstep!(IK::MFPTIntegrationKernel) = eval_driftstep!(IK.sdestep)

function push_Δt!(mfptV, IK::MFPTIntegrationKernel)
    if IK.condition(IK.sdestep.x1)
        push!(mfptV,get_Δt_to_condition(IK))
    else
        push!(mfptV, _Δt(IK.sdestep))
    end
end

function get_Δt_to_condition(IK::MFPTIntegrationKernel{kd,dynT,cT}) where {kd,dynT,cT<:SingleBarrierCondition{i}} where i
    _Δt(IK.sdestep)*(IK.condition.x0 - IK.sdestep.x0[i])/(IK.sdestep.x1[i] - IK.sdestep.x0[i])
end
function get_Δt_to_condition(IK::MFPTIntegrationKernel{kd,dynT,cT}) where {kd,dynT,cT<:DoubleBarrierCondition{i}} where i
    k = IK.condition.condition1(IK.sdestep.x1[i]) ? 1 : 2
    _Δt(IK.sdestep)*(IK.condition.xs[k] - IK.sdestep.x0[i])/(IK.sdestep.x1[i] - IK.sdestep.x0[i])
end

function rescale_discreteintegrator!(IK::MFPTIntegrationKernel{1,dyn}; int_limit_thickness_multiplier = 6, restrict_limit_to_interpolationgrid = true, kwargs...) where dyn <:SDEStep{d,d,m} where {d,m}
    σ = sqrt(_Δt(IK.sdestep)*get_g(IK.sdestep.sde)(d, IK.sdestep.x0,_par(IK.sdestep),_t0(IK.sdestep))^2)
    mn = IK.sdestep.x1[d] - int_limit_thickness_multiplier*σ
    mx = IK.sdestep.x1[d] + int_limit_thickness_multiplier*σ
    if restrict_limit_to_interpolationgrid
        mn = min(IK.T.axes[d][end], max(IK.T.axes[d][1],mn))
        mx = max(IK.T.axes[d][1],min(IK.T.axes[d][end],mx))
    end
    rescale_to_limits!(IK.discreteintegrator, mn, mx)
end


function get_IK_weights!(IK::MFPTIntegrationKernel; kwargs...)
    IK.discreteintegrator(IK, IK.temp.itpM; kwargs...)
end

function (IK::MFPTIntegrationKernel{dk})(vals,x) where dk
    # * if m ≠ 1 and d ≠ k: figure out a rework
    set_oldvals_tozero!(vals, IK)
    update_relevant_states!(IK,x)
    fx = transitionprobability(IK.sdestep,IK.x1)

    if isapprox(fx,zero(fx), atol=1e-8) || IK.condition(IK.x1)
        all_zero!(vals, IK)
    else
        basefun_vals_safe!(IK, allow_extrapolation = true)
        fill_vals!(vals,IK,fx,_idx_it(IK), _val_it(IK))
    end
end

# function basefun_vals_safe!(IK::MFPTIntegrationKernel{dk,sdeT}; kwargs...) where sdeT<:AbstractSDEStep{d} where {dk,d}
#     for (it,ax,x1) in zip(IK.temp.itpVs,IK.T.axes,IK.x1)
#         basefun_vals_safe!(it,ax,x1; allow_extrapolation = true, kwargs..., IK.kwargs...)
#     end
#     nothing
# end

# function update_relevant_states!(IK::MFPTIntegrationKernel{dk,sdeT},x) where sdeT<:SDEStep{d,d,m} where {dk,d,m}
#     @inbounds IK.x1[d] = x
# end
function update_relevant_states!(IK::MFPTIntegrationKernel{dk,sdeT},x) where sdeT<:SDEStep{d,k,m} where {dk,d,k,m}
    for (i,j) in enumerate(k:d)
        @inbounds IK.x1[j] = x[i]
    end
end

function basefun_vals_safe!(IK::MFPTIntegrationKernel{dk,sdeT}; kwargs...) where sdeT<:AbstractSDEStep{d} where {dk,d}
    for (it,ax,x1) in zip(IK.temp.itpVs,IK.T.axes,IK.x1)
        basefun_vals_safe!(it,ax,x1; kwargs..., IK.kwargs...)
    end
    nothing
end

function set_oldvals_tozero!(vals, IK::MFPTIntegrationKernel{kd,dynT,cT,diT,x0T,x1T,TT}) where {kd,dynT,cT,diT,x0T,x1T,TT<:InterpolatedFunction{T,N, itp_type}} where {T,N,itp_type <: SparseInterpolationType}
    @inbounds for idx in _idx_it(IK)
        vals[idx...] = zero(eltype(vals))
    end
end

function all_zero!(vals::AbstractArray{<:T}, IK::MFPTIntegrationKernel{kd,dynT,cT,diT,x0T,x1T,TT}) where {kd,dynT,cT,diT,x0T,x1T,TT<:InterpolatedFunction{T,N, itp_type}} where {T,N,itp_type <: SparseInterpolationType}
end