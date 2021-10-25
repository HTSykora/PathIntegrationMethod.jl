function LinearInterpolation()
    LinearInterpolation(nothing)
end

function Base.getindex(vals::SparseInterpolationBaseVals{1,vT}, idx) where vT
    if idx == vals.idxs[1]
        return vals.val[1]
    elseif idx == vals.idxs[2]
        return vals.val[2]
    end
    return zero(eltype(vT))
end

function linearinterpolation_weights!(vals,x1,x2, x, Δx)
    if Δx isa Nothing
        δ = (x - x1)/(x2-x1);
    else
        δ = (x-x1)/Δx
    end
    vals[1] = one(δ)-δ
    vals[2] = δ
    nothing
end

function basefun_vals!(vals,itp::LinearInterpolation,xs::Vx,x) where {Vx<:AbstractVector{Tx}} where Tx<:Number
    basefun_vals_safe!(vals, itp, xs, x)
    nothing
end

function basefun_vals_safe!(vals::SparseInterpolationBaseVals{ord,vT},itp::LinearInterpolation,xs::Vx,x; kwargs...) where {ord,vT, Vx<:AbstractVector{Tx}} where Tx<:Number
    do_interpolation, zero_extrapolation, i = find_idx(xs, x,; kwargs...)
    # for j in eachindex(vals.val)
    #     vals.val[j] = zero(eltype(vT))
    # end

    if do_interpolation
        vals.idxs[1] = i
        vals.idxs[2] = i+1
        linearinterpolation_weights!(vals.val, xs[i], xs[i+1], x, itp.Δ)
    elseif zero_extrapolation
        vals.val .= zero(eltype(vT))
    elseif i==1
        vals.val[1] = one(eltype(vT))
        vals.val[2] = zero(eltype(vT))
        vals.idxs[1] = 1
        vals.idxs[2] = 2
    else
        vals.val[1] = zero(eltype(vT))
        vals.val[2] = one(eltype(vT))
        vals.idxs[1] = i-1
        vals.idxs[2] = i
    end
        
    nothing
    # return vals
end