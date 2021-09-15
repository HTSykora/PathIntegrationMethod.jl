# Constructor
function ChebyshevInterpolation(N::Integer;)
    ChebyshevInterpolation{N-1}()
end

##  Functions for the interpolation and integration
ξ2x(ξ, a, b) = (ξ + 1)*((b - a)/2) + a # local to global
x2ξ(x, a, b) = 2*(x - a)/(b - a) - 1 # global to local
x2θ(x, a, b) = acos(x2ξ(x, a, b)) # global to angular

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

function _b(T::DataType,n,j)
    if j == n÷2
        one(T)
    else
        T(2)
    end
end
"""
    clenshawcurtisweights([T,],n)
Create an array of `n` weights corresponding to the `n` chebyshev nodes in [-1,1] with type T or Float64
"""
clenshawcurtisweights(n1::Integer) = clenshawcurtisweights(Float64,n1::Integer)
function clenshawcurtisweights(T,n1::Integer)
    n = n1-1
    _n = 1/n
    w = zeros(T,n1);
    w[1] = (1-sum(_b(T,n,j)/(4j^2-1) for j in 1:n÷2))*_n
    w[end] = w[1]#_w0(n1)#(1-sum(_b(Float64,n,j)/(4j^2-1) for j in 1:n÷2))*_n
    for k in 1:n-1
        ϑk = k*_n
        w[k+1] = 2*(1-sum(_b(T,n,j)/(4j^2-1)*cospi(2j*ϑk) for j in 1:n÷2))*_n
    end
    w
end
clenshawcurtisweights(xa::Number,xb::Number,n1::Integer) = clenshawcurtisweights(Float64,xa,xb,n1)
function clenshawcurtisweights(T,xa::Number,xb::Number,n1::Integer)
    0.5*(xb-xa) .* clenshawcurtisweights(T,n1)
end

# function (itp::ChebyshevInterpolation{N})(p::Vp,xs::Vx,x::xT, i::Integer) where {N,Vp<:AbstractVector{Tp},Vx<:AbstractVector{Tx},xT<:Number} where Tp<:Number where Tx<:Number
#     _1 = iseven(N) ? 1 : -1
#     s1 = 0.5*(p[1]/(x-xs[1]) +  _1 * p[N+1]/(x-xs[N+1]))
#     s2 = 0.5*(1/(x-xs[1]) +  _1 /(x-xs[N+1]))
#     _1 = 1
#     for j in 2:N
#         _1 *= -1
#         _1px = 1/(x-xs[j])
#         s1 += _1*p[j]*_1px
#         s2 += _1*_1px
#     end
#     return s1/s2
# end

function basefun_vals!(itp::ChebyshevInterpolation{N},vals,xs::Vx,x) where {N,Vx<:AbstractVector{Tx}} where Tx<:Number
    # i = 0... N
    _1 = iseven(N) ? 1 : -1
    vals[1] = 0.5/(x-xs[1])
    vals[end] =  0.5* _1 /(x-xs[N+1])
    s2 = 0.5*(1/(x-xs[1]) +  _1 /(x-xs[N+1]))
    _1 = 1
    for j in 2:N
        _1 *= -1
        _1px = 1/(x-xs[j])
        vals[j] = _1*_1px
        s2 += _1*_1px
    end
    vals .= vals./s2
    return vals
end

function basefun_vals_safe!(itp::ChebyshevInterpolation{N},vals,xs::Vx,x; allow_extrapolation = false, kwargs...) where {N,Vx<:AbstractVector{Tx}} where Tx<:Number
    needsinterpolation, i = find_idx(xs, x, allow_extrapolation = allow_extrapolation; kwargs...)
    if needsinterpolation
        basefun_vals!(itp,vals,xs,x)
    else
        vals .= zero(eltype(vals))
        vals[i] = one(eltype(vals))
    end
    return vals
end