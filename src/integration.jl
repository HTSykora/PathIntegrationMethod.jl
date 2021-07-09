function _integrate(p::Vp,xs::NTuple{1,xsT}) where {Vp<:AbstractArray{Tp,1},xsT<:Axis{itpT}} where {Tp<:Number, itpT<:LinearInterpolation{true}}
    _integrate(p,xs...)
end
function _integrate(p::Vp,xs::NTuple{1,xsT}) where {Vp<:AbstractArray{Tp,1},xsT<:Axis} where {Tp<:Number}
    _integrate(p,xs...)
end

function _integrate(f::Function,p::Vp,xs::NTuple{1,xsT}) where {Vp<:AbstractVector{Tp},xsT<:Axis} where Tp<:Number 
    _integrate(f,p,xs...)
end


function _integrate(p::Vp,x::Axis{itpT}) where {Vp<:AbstractVector{Tp},itpT<:LinearInterpolation{true}} where Tp<:Number 
    x.itp.Δ*(sum(p[i] for i in 2:(length(p)-1)) + 0.5*(p[1]+p[end]))
end
function _integrate(f::Function, p::Vp,x::Axis{itpT}) where {Vp<:AbstractVector{Tp},itpT<:LinearInterpolation{true}} where Tp<:Number 
    x.itp.Δ*(sum(p[i]*f(x[i]) for i in 2:(length(p)-1)) + 0.5*(p[1]*f(x[1])+p[end]*f(x[end])))
end

function _integrate(p::Vp,xs::NTuple{N,xsT}) where {Vp<:AbstractArray{Tp,N},xsT<:Axis{itpT}} where {N,Tp<:Number,itpT<:LinearInterpolation{true}}
    sp = size(p);
    last(xs).itp.Δ* sum(last(xs).itp.wts[i]*_integrate(view(p,(Colon() for _ in 1:N-1)...,i),xs[1:N-1]) for i in 1:sp[end])
    # ΔA = prod(x.itp.Δ for x in xs)
    # corners = Base.Iterators.product(prepend_1.(sp)...)
    # S = sum(p[corner...] for corner in corners) * 0.5^length(sp);
    # *(sum(p[i] for i in 2:(length(p)-1)) + 0.5*(p[1]+p[end]))
end


function _integrate(p::Vp,x::Axis) where {Vp<:AbstractVector{Tp}} where Tp<:Number 
    sum(x.itp.wts[i]*_p for (i,_p) in enumerate(p))
end
function _integrate(f::Function, p::Vp,x::Axis) where {Vp<:AbstractVector{Tp}} where Tp<:Number 
    sum(x.itp.wts[i]*_p*f(x[i]) for (i,_p) in enumerate(p))
end