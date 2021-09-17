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

get_axistempval(axis,i) = axis.temp[i]
get_axistempval(axes,idx,i) = axes[i].temp[idx[i]]
reduce_axestempprod(ax1) = one(eltype(ax1.temp))
function reduce_axestempprod(axi1,axi_tail...)
    val = get_axistempval(axi1...)
    val * reduce_axestempprod(axi_tail...)
end

# Multivariate interopolation
function interpolate(p::MX,axes,x::Vararg{Any,N}; idx_it = Base.Iterators.product(eachindex.(axes)...)) where MX<:AbstractArray{T,N} where {T,N} #xs <: NTuple{<:gridAxis}
    val = zero(eltype(p))
    # @inbounds for (i,_x) in enumerate(x)
    #     basefun_vals_safe!(axes[i],_x)
    # end
    #
    basefun_vals_safe!.(axes,x)
    # basefun_vals_safe!(axes[1],x[1])
    # basefun_vals_safe!(axes[2],x[2])
    # basefun_vals_safe!(axes[3],x[3])
    
    @inbounds for idx in idx_it
        # _temp = one(val)
        # @inbounds for i in 1:N
        #     # axes[i].temp[idx[i]]
        #     _temp *= get_axistempval(axes,idx,i)
        # end
        # _temp *= axes[1].temp[idx[1]]
        # _temp *= axes[2].temp[idx[2]]
        # _temp *= axes[3].temp[idx[3]]
        _temp = mapreduce(get_axistempval, *, axes, idx)
        # _temp = prod(a.temp[idx[i]] for (i,a) in enumerate(axes))
        # _temp = reduce_axestempprod(zip(axes,idx)...)
        val += p[idx...] * _temp
        # val += p[idx...] * prod(axes[i].temp[_i] for (i,_i) in enumerate(idx))
        # val += p[idx...] * prod(axes[i].temp[_i] for (i,_i) in enumerate(idx))
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
