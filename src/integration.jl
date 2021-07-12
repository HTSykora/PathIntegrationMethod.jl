# Integration kernel for fast computation of the transition matrix
IntegrationKernel(idx₀::iT,idx₁::iT,xs::xT,f::fT,TT::ttT) where {iT,xT,fT,ttT} = 
    IntegrationKernel{iT,iT,xT,fT,ttT,Nothing}(idx₀,idx₁,xs,f,TT,nothing)

function (b::IntegrationKernel)(vals,x₀)
    basefun_vals!(b.xs.itp,vals,b.xs,x₀)
    vals .*= _tp(b.TT.sde,b.TT.pdgrid.xs[1][b.idx₁...],b.TT.Δt, x₀,zero(b.TT.Δt), method = b.TT.method)
    vals
end
function (b::IntegrationKernel)(x₀)
    basefun_vals(b.xs.itp,b.xs,x₀) .* _tp(b.TT.sde,b.TT.pdgrid.xs[1][b.idx₁...],b.TT.Δt, x₀,zero(b.TT.Δt), method = b.TT.method)
end

# Integration over the axis

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
    x.wts.Δ*(sum(p[i] for i in 2:(length(p)-1)) + 0.5*(p[1]+p[end]))
end
function _integrate(f::Function, p::Vp,x::Axis{itpT}) where {Vp<:AbstractVector{Tp},itpT<:LinearInterpolation{true}} where Tp<:Number 
    x.wts.Δ*(sum(p[i]*f(x[i]) for i in 2:(length(p)-1)) + 0.5*(p[1]*f(x[1])+p[end]*f(x[end])))
end

function _integrate(p::Vp,x::Axis) where {Vp<:AbstractVector{Tp}} where Tp<:Number 
    sum(x.wts[i]*_p for (i,_p) in enumerate(p))
end
function _integrate(f::Function, p::Vp,x::Axis) where {Vp<:AbstractVector{Tp}} where Tp<:Number 
    sum(x.wts[i]*_p*f(x[i]) for (i,_p) in enumerate(p))
end


# TODO: function integration over the axis
function _integrate(p::Vp,xs::NTuple{N,xsT}) where {Vp<:AbstractArray{Tp,N},xsT<:Axis{itpT}} where {N,Tp<:Number,itpT<:LinearInterpolation{true}}
    sp = size(p);
    last(xs).wts.Δ* sum(last(xs).wts[i]*_integrate(view(p,(Colon() for _ in 1:N-1)...,i),xs[1:N-1]) for i in 1:sp[end])
end
function _integrate(p::Vp,xs::Tuple{xsT}) where {Vp<:AbstractArray{Tp,N},xsT<:Axis{itpT}} where {N,Tp<:Number}
    sp = size(p,N);
    sum(last(xs).wts[i]*_integrate(view(p,(Colon() for _ in 1:N-1)...,i),xs[1:N-1]) for i in 1:sp)
end


