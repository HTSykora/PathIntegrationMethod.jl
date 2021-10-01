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
function clenshawcurtisweights(T::DataType,n1::Integer)
    # N1 = 10001; b = 1; a = -1;
    # bma=(b-a);
    N=n1-1; 
    Ns = 2:2:N;
    _w= zeros(T,2*(n1-1));
    w = zeros(T,n1);

    _w[1] = T(2);
    for (i,ic) in enumerate(3:2:n1)
        _w[ic] = T(2)/(one(T)-T(Ns[i])^2)
    end
    for ic in 3:2:n1-1
        _w[end-ic+2] = _w[ic]
    end
    _iw = real.(ifft(_w))
    
    w[1] = _iw[1]
    w[end] = _iw[n1]
    for i in 2:length(w)-1
        w[i] = T(2)*_iw[i]
    end
    w
end
function clenshawcurtisweights(T::DataType,xa,xb,n1::Integer)
    bma2 = T(0.5)*(T(xb)-T(xa))
    w = clenshawcurtisweights(T,n1)
    w .*= bma2
    w
end
function clenshawcurtisweights(n1::Integer) where {N}
    clenshawcurtisweights(Float64,n1)
end
function clenshawcurtisweights(xa::Ta,xb::Tb,n1::Integer) where {N,Ta<:Number,Tb<:Number}
    clenshawcurtisweights(Float64,xa,xb,n1)
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