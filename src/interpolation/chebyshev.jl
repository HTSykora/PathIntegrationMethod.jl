
ξ2x(ξ, a, b) = (ξ + 1)*((b - a)/2) + a
x2ξ(x, a, b) = 2*(x - a)/(b - a) - 1
x2θ(x, a, b) = acos(x2ξ(x, a, b))

"""
    chebygrid(n)
Create an array of `n` chebyshev nodes in [-1,1]
"""
function chebygrid(n::Integer)
    cos.(π*(n-1:-1:0)/(n-1))
end

"""
    chebygrid(xa, xb, n)
Create an array of `n` chebyshev nodes in [`xa`,`xb`]
"""
function chebygrid(xa, xb, n::Integer)
    ξ2x.(chebygrid(n), xa, xb)
end

function (itp::ChebyshevInterpolation{N,1})(p::Vp,xs::Vx,x::xT) where {N,Vp<:AbstractVector{Tp},Vx<:AbstractVector{Tx},xT<:Number} where Tp<:Number where Tx<:Number
    if !(x<xs[1] || x>xs[end])
        if x isa Union{Rational,Integer}
            @inbounds i = searchsortedlast(xs, x)  
        else
            @inbounds i = searchsortedlast(xs, x - 10eps(xT))
        end
        if i==0  # it can happen if t ≈ t0
            return p[1]
        elseif abs2(xs[i+1] - x) < eps()
            return p[i+1]
        end
    else
        _1 = iseven(N) ? 1 : -1
        s1 = 0.5*(p[1]/(x-xs[1]) +  _1 * p[N+1]/(x-xs[N+1]))
        s2 = 0.5*(1/(x-xs[1]) +  _1 /(x-xs[N+1]))
        _1 = 1
        for j in 2:N
            _1 *= -1
            _1px = 1/(x-xs[j])
            s1 += _1*p[j]*_1px
            s2 += _1*_1px
        end
        return s1/s2
    end
end