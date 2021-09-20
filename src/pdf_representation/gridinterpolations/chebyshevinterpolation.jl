# Constructor
function ChebyshevInterpolation(N::Integer;)
    _N = N - 1
    _1 = iseven(_N) ? 1 : -1
    ChebyshevInterpolation(_N,_1)
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
    chebygrid(Float64, n)
end
function chebygrid(T, n::Integer)
    T.(cos.(π*(n-1:-1:0)/(n-1)))
end

"""
    chebygrid(xa, xb, n)
Create an array of `n` chebyshev nodes in [`xa`,`xb`]
"""
function chebygrid(T, xa, xb, n::Integer)
    ξ2x.(chebygrid(T, n), T(xa), T(xb))
end
function chebygrid(xa, xb, n::Integer)
    chebygrid(Float64, xa, xb, n)
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
function clenshawcurtisweights(T, n1::Integer)
    n = n1-1
    _n = 1/n
    w = zeros(T,n1);
    w[1] = (1-sum(_b(T,n,j)/(4j^2-1) for j in 1:n÷2))*_n
    w[end] = w[1]#_w0(n1)#(1-sum(_b(Float64,n,j)/(4j^2-1) for j in 1:n÷2))*_n
    @inbounds for k in 1:n-1
        ϑk = k*_n
        w[k+1] = 2*(1-sum(_b(T,n,j)/(4j^2-1)*cospi(2j*ϑk) for j in 1:n÷2))*_n
    end
    w
end
clenshawcurtisweights(xa::Number,xb::Number,n1::Integer) = clenshawcurtisweights(Float64,xa,xb,n1)
function clenshawcurtisweights(T,xa::Number,xb::Number,n1::Integer)
    T(0.5)*(T(xb)-T(xa)) .* clenshawcurtisweights(T,n1)
end

function basefun_vals!(vals,itp::ChebyshevInterpolation,xs::Vx,x) where {Vx<:AbstractVector{Tx}} where Tx<:Number
    # i = 0... N
    N = itp.N
    _1 = itp._1 # iseven(N) ? 1 : -1
    vals[1] = 0.5/(x-xs[1])
    vals[end] =  0.5* _1 /(x-xs[N+1])
    s2 = 0.5*(1/(x-xs[1]) +  _1 /(x-xs[N+1]))
    _1 = 1
    @inbounds for j in 2:N
        _1 *= -1
        _1px = 1/(x-xs[j])
        vals[j] = _1*_1px
        s2 += _1*_1px
    end
    @inbounds for i in eachindex(vals)
        vals[i] = vals[i]/s2
    end
    # vals .= vals ./ s2
    nothing
    # return vals
end

function basefun_vals_safe!(vals,itp::ChebyshevInterpolation,xs::Vx,x; allow_extrapolation = false, kwargs...) where {Vx<:AbstractVector{Tx}} where Tx<:Number
    needsinterpolation, i = find_idx(xs, x, allow_extrapolation = allow_extrapolation; kwargs...)
    if needsinterpolation
        basefun_vals!(vals,itp,xs,x)
    else
        for j in eachindex(vals)
            vals[j] = zero(eltype(vals))
        end
        vals[i] = one(eltype(vals))
    end
    nothing
    # return vals
end