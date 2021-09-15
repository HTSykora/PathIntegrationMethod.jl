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


function find_idx(xs::Vx,x::xT; allow_extrapolation::Bool = false, zero_extrapolation = false, kwargs...) where {N,Vx<:AbstractVector{Tx},xT<:Number} where Tp<:Number where Tx<:Number
    # this function is also used in basefun_vals_safe! (Chebyshev)
    if !(x<xs[1] || x>xs[end])
        if x isa Union{Rational,Integer}
            @inbounds i = searchsortedlast(xs, x)  
        else
            @inbounds i = searchsortedlast(xs, x - 10eps(xT))
        end
        if i==0  # it can happen if t ≈ t0
            return false, one(typeof(i))
        elseif isapprox(xs[i+1], x, atol = 100eps(xT))
            return false, i+one(typeof(i))
        else
            return true, i
        end
    else
        i = x<xs[1] ? 1 : length(xs);
        return allow_extrapolation, i
    end
end

# Multivariate interopolation
function interpolate(p::MX,axes,x::Vararg{Any,N}) where MX<:AbstractArray{T,N} where {T,N} #xs <: NTuple{<:gridAxis}
    val = zero(eltype(p))
    for (i,_x) in enumerate(x)
        basefun_vals_safe!(axes[i],_x)
    end
    # c_idxs = CartesianIndices(size(p))
    idx_it = Base.Iterators.product(eachindex.(axes)...)
    for idx in idx_it
        val += p[idx...] * prod(axes[i].temp[_i] for (i,_i) in enumerate(idx))
    end
    
    val
end

function basefun_vals_safe!(vals,axis::GridAxis,x)
    basefun_vals_safe!(vals,axis.itp,axis.xs,x)
end
function basefun_vals_safe!(axis::GridAxis,x)
    basefun_vals_safe!(axis.temp,axis.xs,x)
end
function basefun_vals_safe(axis,x) where Vx<:AbstractVector{Tx} where Tx<:Number
    # i = 0... N
    vals = zeros(length(axis));
    basefun_vals_safe!(vals,axis,x)
    return vals
end
