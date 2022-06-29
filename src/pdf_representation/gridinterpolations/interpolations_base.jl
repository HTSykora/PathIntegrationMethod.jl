################################
## Generic interpolation functions

# 1D interpolation
# function (itp::AbstractInterpolationType)(p::Vp,xs::Vx,x::xT) where {N,Vp<:AbstractVector{Tp},Vx<:AbstractVector{Tx},xT<:Number} where Tp<:Number where Tx<:Number
#     do_interpolation, zero_extrapolation, i = find_idx(xs, x)
#     if do_interpolation
#         return itp(p,xs,x,i)
#     else
#         return p[i]
#     end
# end


function find_idx(xs::Vx,x::xT; allow_extrapolation::Bool = false, zero_extrapolation::Bool = true, kwargs...) where {Vx<:AbstractVector{Tx},xT<:Number} where Tx<:Number
    # this function is also used in basefun_vals_safe! (Chebyshev)
    if !(x<xs[1] || x>xs[end])
        if x isa Union{Rational,Integer}
            @inbounds i = searchsortedlast(xs, x)  
        else
            @inbounds i = searchsortedlast(xs, x - 10eps(xT))
        end
        if i==zero(i)  # it can happen if t ≈ t0
            do_interpolation = false
            i = one(i)
        elseif i == length(xs)
            do_interpolation = false
        elseif isapprox(xs[i+1], x, atol = 100eps(xT))
            do_interpolation = false
            i = i + one(i)
        else
            do_interpolation = true
        end
        _ze = false
    else
        do_interpolation = allow_extrapolation
        _ze = zero_extrapolation && !(allow_extrapolation)
        i = x<xs[1] ? 1 : length(xs);
    end
    return do_interpolation, _ze, i
end

# Multivariate interopolation
get_tempval(axis::AxisGrid,i) = axis.temp[i]
function interpolate(p::MX,axes, x::Vararg{Any,N}; idx_it = BI_product(_eachindex.(axes)...), val_it = BI_product(_gettempvals.(axes)...), kwargs...) where MX<:AbstractArray{T,N} where {T,N} #xs <: NTuple{<:gridAxis}
    for (axis, _x) in zip(axes,x)
        basefun_vals_safe!(axis,_x; kwargs...)
    end
    val = zero(eltype(p))
    @inbounds for (idx, _val) in zip(idx_it, val_it)
        val += p[idx...] * prod(_val)#reduce_tempprod(zip(axes,idx)...)
    end
    
    val
end

function basefun_vals_safe!(vals,axis::AxisGrid,x; kwargs...)
    basefun_vals_safe!(vals,axis.itp,axis.xs,x; kwargs...)
    nothing
end
function basefun_vals_safe!(axis::AxisGrid,x; kwargs...)
    basefun_vals_safe!(axis.temp, axis,x; kwargs...)
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
    idxs = @MVector ones(lT,order +1)
    val = @MVector zeros(T,length(idxs))
    SparseInterpolationBaseVals{order,typeof(val),typeof(idxs),lT}(val,idxs,l)
end
Base.size(vals::SparseInterpolationBaseVals) = (vals.l,)
Base.size(vals::SparseInterpolationBaseVals,i) = i==1 ? vals.l : one(vals.l)
Base.length(vals::SparseInterpolationBaseVals) = vals.l
_val(vals::SparseInterpolationBaseVals) = vals.val
_eachindex(vals::SparseInterpolationBaseVals) where itpT<:LinearInterpolation = vals.idxs

Base.similar(lbfv::SparseInterpolationBaseVals) = SparseInterpolationBaseVals(order, similar(lbfv.val),similar(lbfv.idxs),lbfv.l)
Base.zero(lbfv::SparseInterpolationBaseVals{order}) where order = SparseInterpolationBaseVals(order, zero(lbfv.val),one.(lbfv.idxs),lbfv.l)
