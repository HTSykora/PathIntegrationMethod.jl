abstract type AbstractInterpolationType{TF} end

struct ChebyshevInterpolation{N,Ttmp} <: AbstractInterpolationType{false}
    tmp::Ttmp
end
ChebyshevInterpolation(N) = ChebyshevInterpolation{N,Nothing}(nothing)

struct LinearInterpolation{TF,Ttmp,tΔ} <: AbstractInterpolationType{TF}
    tmp::Ttmp
    Δ::tΔ
end

LinearInterpolation(lvl; Q_equidistant = false) = LinearInterpolation{Q_equidistant,Nothing,Nothing}(nothing,nothing)
function LinearInterpolation(lvl, Δx; Q_equidistant = false) 
    _Δx = Q_equidistant ? Δx : nothing
    LinearInterpolation{Q_equidistant,Nothing,typeof(_Δx)}(nothing,_Δx)
end



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
            return true, i+1, zero(Tp)
        end
    else
        return true, 0, zero(Tp)
    end
end
