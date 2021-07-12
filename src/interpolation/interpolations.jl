# Trapezoidal weights for linear interpolation
Base.getindex(t::TrapezoidalWeights{T},idx::Integer) where T = (idx == 1 || idx == t.l) ? T(0.5) : T(1.)
Base.size(t::TrapezoidalWeights) = (t.l,)

################################
## Generic functions
function (itp::AbstractInterpolationType)(p::Vp,xs::Vx,x::xT) where {N,Vp<:AbstractVector{Tp},Vx<:AbstractVector{Tx},xT<:Number} where Tp<:Number where Tx<:Number
    needsinterpolation, i, val = find_idx(p, xs, x)
    if needsinterpolation
        return itp(p,xs,x,i)
    else
        return val 
    end
end

function  find_idx(p::Vp,xs::Vx,x::xT) where {N,Vp<:AbstractVector{Tp},Vx<:AbstractVector{Tx},xT<:Number} where Tp<:Number where Tx<:Number
    if !(x<xs[1] || x>xs[end])
        if x isa Union{Rational,Integer}
            @inbounds i = searchsortedlast(xs, x)  
        else
            @inbounds i = searchsortedlast(xs, x - 10eps(xT))
        end
        if i==0  # it can happen if t ≈ t0
            return false, 1, p[1]
        elseif abs2(xs[i+1] - x) < eps()
            return false, i+1, p[i+1]
        else
            return true, i, zero(Tp)
        end
    else
        return true, 0, zero(Tp)
    end
end