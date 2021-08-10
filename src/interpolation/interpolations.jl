# Trapezoidal weights for linear interpolation
Base.getindex(t::TrapezoidalWeights{T},idx::Integer) where T = (idx == 1 || idx == t.l) ? T(0.5) : T(1.)
Base.size(t::TrapezoidalWeights) = (t.l,)


################################
## Generic interpolation functions

# 1D interpolation
function (itp::AbstractInterpolationType)(p::Vp,xs::Vx,x::xT) where {N,Vp<:AbstractVector{Tp},Vx<:AbstractVector{Tx},xT<:Number} where Tp<:Number where Tx<:Number
    needsinterpolation, i = find_idx(xs, x)
    if needsinterpolation
        return itp(p,xs,x,i)
    else
        return p[i]
    end
end


function find_idx(xs::Vx,x::xT; allow_extrapolation = false, zero_extrapolation = false, kwargs...) where {N,Vx<:AbstractVector{Tx},xT<:Number} where Tp<:Number where Tx<:Number
    # this function is also used in basefun_vals_safe! (Chebyshev)
    if !(x<xs[1] || x>xs[end])
        if x isa Union{Rational,Integer}
            @inbounds i = searchsortedlast(xs, x)  
        else
            @inbounds i = searchsortedlast(xs, x - 10eps(xT))
        end
        if i==0  # it can happen if t ≈ t0
            return false, 1
        elseif abs2(xs[i+1] - x) < eps()
            return false, i+1
        else
            return true, i
        end
    else
        return allow_extrapolation, x<xs[1] ? 1 : length(xs)
    end
end

# Multivariate interopolation
function interpolate_MV(p::MX,xs,x...) where MX<:AbstractArray{T,N} where {T,N} where AX<:Axis
    val = zero(eltype(p))
    for (i,_x) in enumerate(x)
        basefun_vals_safe!(xs[i].itp.tmp,xs[i],_x)
    end
    # c_idxs = CartesianIndices(size(p))
    idx_it = Base.Iterators.product(eachindex.(xs)...)
    for idx in idx_it
        val += p[idx...] * prod(xs[i].itp.tmp[_i] for (i,_i) in enumerate(idx))
    end
    
    val
end


function basefun_vals!(vals,ax::Axis,x)
    basefun_vals!(ax.itp,vals,ax,x)
end

function basefun_vals(itp,xs::Vx,x) where Vx<:AbstractVector{Tx} where Tx<:Number
    # i = 0... N
    vals = zeros(N+1);
    basefun_vals!(itp,vals,xs,x)
    return vals
end

function basefun_vals(ax::Axis,x)
    basefun_vals(ax.itp,ax,x)
end


function basefun_vals_safe!(vals,ax::Axis,x)
    basefun_vals_safe!(ax.itp,vals,ax,x)
end

function basefun_vals_safe(itp,xs::Vx,x) where Vx<:AbstractVector{Tx} where Tx<:Number
    # i = 0... N
    vals = zeros(length(xs));
    basefun_vals_safe!(itp,vals,xs,x)
    return vals
end
function basefun_vals_safe(ax::Axis,x)
    basefun_vals_safe(ax.itp,ax,x)
end
