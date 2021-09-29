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
    val = zeros(T,2)
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
Base.size(vals::LinearBaseFunVals,i) = i==1 ? vals.l : one(vals.l)
Base.length(vals::LinearBaseFunVals) = vals.l
_val(vals::LinearBaseFunVals) = vals.val
_eachindex(vals::LinearBaseFunVals) where itpT<:LinearInterpolation = vals.idxs

_eachindex(axis::GridAxis{itpT}) where itpT<:LinearInterpolation = _eachindex(axis.temp)
_gettempvals(axis::GridAxis{itpT}) where itpT<:LinearInterpolation = axis.temp.val
Base.similar(lbfv::LinearBaseFunVals) = LinearBaseFunVals(similar(lbfv.val),similar(lbfv.idxs),lbfv.l)
Base.zero(lbfv::LinearBaseFunVals) = LinearBaseFunVals(zero(lbfv.val),one.(lbfv.idxs),lbfv.l)

function duplicate(a::GridAxis{itpT,wT,xT,xeT,tmpT}) where {itpT<:SparseInterpolationType,wT,xT,xeT,tmpT<:LinearBaseFunVals} 
    GridAxis{itpT,wT,xT,xeT,tmpT}(deepcopy(a.itp), deepcopy(a.xs), deepcopy(a.wts), zero(a.temp))
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

function basefun_vals_safe!(vals::LinearBaseFunVals{vT},itp::LinearInterpolation,xs::Vx,x; allow_extrapolation = false, kwargs...) where {vT, Vx<:AbstractVector{Tx}} where Tx<:Number
    needsinterpolation, i = find_idx(xs, x, allow_extrapolation = allow_extrapolation; kwargs...)
    for j in eachindex(vals.val)
        vals.val[j] = zero(eltype(vT))
    end

    if needsinterpolation
        vals.idxs[1] = i
        vals.idxs[2] = i+1
        linearinterpolation_weights!(vals.val, xs[i], xs[i+1], x, itp.Δ)
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