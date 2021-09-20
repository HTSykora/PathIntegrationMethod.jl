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
        if i==0  # it can happen if t ≈ t0
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
function interpolate(p::MX,axes, x::Vararg{Any,N}; idx_it = Base.Iterators.product(eachindex.(axes)...)) where MX<:AbstractArray{T,N} where {T,N} #xs <: NTuple{<:gridAxis}
    basefun_vals_safe!.(axes,x)
    val = zero(eltype(p))
    @inbounds for idx in idx_it
        val += p[idx...] * reduce_tempprod(zip(axes,idx)...)
    end
    
    val
end

function basefun_vals_safe!(vals,axis::GridAxis,x)
    basefun_vals_safe!(vals,axis.itp,axis.xs,x)
end
function basefun_vals_safe!(axis::GridAxis,x)
    basefun_vals_safe!(axis.temp, axis,x)
end
function basefun_vals_safe(axis,x) where Vx<:AbstractVector{Tx} where Tx<:Number
    # i = 0... N
    vals = zeros(length(axis));
    basefun_vals_safe!(vals,axis,x)
    return vals
end
