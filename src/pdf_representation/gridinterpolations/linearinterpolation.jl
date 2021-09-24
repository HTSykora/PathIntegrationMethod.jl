function LinearInterpolation()
    LinearInterpolation(nothing)
end

function Base.getindex(tw::TrapezoidalWeights, idx) 
    if idx == 1 || idx == tw.l
        return eltype(tw.Δ)(0.5)*tw.Δ
    else
        return tw.Δ
    end
end
Base.size(tw::TrapezoidalWeights) = (tw.l,)

function LinearBaseFunVals(T::DataType, l)
    idxs = [one(l), one(l)]
    val = Vector{T}(undef,2)
    LinearBaseFunVals(val,idxs,l)
end
function Base.getindex(vals::LinearBaseFunVals{vT}, idx) where vT
    if idx == vals.idxs[1]
        return vals.val[1]
    elseif idx == vals.idxs[2]
        return vals.val[2]
    end
    return zero(eltype(vT))
end
Base.size(vals::LinearBaseFunVals) = (vals.l,)
Base.eachindex(axis::GridAxis{itpT}) where itpT<:LinearInterpolation = eachindex(axis.temp)
Base.eachindex(vals::LinearBaseFunVals) = vals.idxs

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

function basefun_vals_safe!(vals::LinearBaseFunVals{vT},itp::LinearInterpolation,xs::Vx,x; allow_extrapolation = false, kwargs...) where {vT, Vx<:AbstractVector{Tx}} where Tx<:Number
    needsinterpolation, i = find_idx(xs, x, allow_extrapolation = allow_extrapolation; kwargs...)
    for j in eachindex(vals.val)
        vals.val[j] = zero(eltype(vT))
    end

    if needsinterpolation
        vals.idxs[1] = i
        vals.idxs[2] = i+1
        linearinterpolation_weights!(vals.val, xs[i], xs[i+1], x, itp.Δ)
    else
        vals.val[1] = one(eltype(vT))
        vals.val[2] = zero(eltype(vT))
        vals.idxs[1] = i
        vals.idxs[2] = i+1
    end
        
    nothing
    # return vals
end