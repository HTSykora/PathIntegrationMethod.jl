################################
## Generic interpolation functions

# 1D interpolation
# function (itp::AbstractInterpolationType)(p::Vp,xs::Vx,x::xT) where {N,Vp<:AbstractVector{Tp},Vx<:AbstractVector{Tx},xT<:Number} where Tp<:Number where Tx<:Number
#     needsinterpolation, i = find_idx(xs, x)
#     if needsinterpolation
#         return itp(p,xs,x,i)
#     else
#         return p[i]
#     end
# end


function find_idx(xs::Vx,x::xT; allow_extrapolation::Bool = false, zero_extrapolation = false, kwargs...) where {Vx<:AbstractVector{Tx},xT<:Number} where Tx<:Number
    # this function is also used in basefun_vals_safe! (Chebyshev)
    if !(x<xs[1] || x>xs[end])
        if x isa Union{Rational,Integer}
            @inbounds i = searchsortedlast(xs, x)  
        else
            @inbounds i = searchsortedlast(xs, x - 10eps(xT))
        end
        if i==zero(i)  # it can happen if t ≈ t0
            needsinterpolation = false
            i = one(i)
        elseif isapprox(xs[i+1], x, atol = 100eps(xT))
            needsinterpolation = false
            i = i + one(i)
        else
            needsinterpolation = true
        end
    else
        needsinterpolation = allow_extrapolation
        i = x<xs[1] ? 1 : length(xs);
    end
    return needsinterpolation, i
end


# Multivariate interopolation
get_tempval(axis::GridAxis,i) = axis.temp[i]
function interpolate(p::MX,axes, x::Vararg{Any,N}; idx_it = BI_product(_eachindex.(axes)...), val_it = BI_product(_gettempvals.(axes)...)) where MX<:AbstractArray{T,N} where {T,N} #xs <: NTuple{<:gridAxis}
    for (axis, _x) in zip(axes,x)
        basefun_vals_safe!(axis,_x)
    end
    val = zero(eltype(p))
    @inbounds for (idx, _val) in zip(idx_it, val_it)
        val += p[idx...] * prod(_val)#reduce_tempprod(zip(axes,idx)...)
    end
    
    val
end

function basefun_vals_safe!(vals,axis::GridAxis,x)
    basefun_vals_safe!(vals,axis.itp,axis.xs,x)
    nothing
end
function basefun_vals_safe!(axis::GridAxis,x)
    basefun_vals_safe!(axis.temp, axis,x)
    nothing
end
function basefun_vals_safe(axis,x) where Vx<:AbstractVector{Tx} where Tx<:Number
    # i = 0... N
    vals = zeros(length(axis));
    basefun_vals_safe!(vals,axis,x)
    return vals
end


# Sparse interpolation
# Getindex: define at sparse interpolation
function SparseInterpolationBaseVals(order, val::valT, idxs::idxsT, l::lT) where {valT, idxsT, lT}
    SparseInterpolationBaseVals{order,valT,idxsT,lT}(val,idxs,l)
end
function SparseInterpolationBaseVals(T::DataType, l::lT, order = 1) where lT<:Integer
    idxs = ones(lT,order +1)
    val = zeros(T,length(idxs))
    SparseInterpolationBaseVals{order,typeof(val),typeof(idxs),lT}(val,idxs,l)
end
Base.size(vals::SparseInterpolationBaseVals) = (vals.l,)
Base.size(vals::SparseInterpolationBaseVals,i) = i==1 ? vals.l : one(vals.l)
Base.length(vals::SparseInterpolationBaseVals) = vals.l
_val(vals::SparseInterpolationBaseVals) = vals.val
_eachindex(vals::SparseInterpolationBaseVals) where itpT<:LinearInterpolation = vals.idxs

Base.similar(lbfv::SparseInterpolationBaseVals) = SparseInterpolationBaseVals(order, similar(lbfv.val),similar(lbfv.idxs),lbfv.l)
Base.zero(lbfv::SparseInterpolationBaseVals{order}) where order = SparseInterpolationBaseVals(order, zero(lbfv.val),one.(lbfv.idxs),lbfv.l)
