# Integration kernel for fast computation of the transition matrix
function IntegrationKernel(sde::sdeT,f::fT,xs::xT,idx₀::iT0,idx₁::iT1,pdgrid::pdT,t₀::tT,t₁::tT,method::methodT) where {sdeT, fT, xT, iT0, iT1, pdT, tT, methodT}
    IntegrationKernel{sdeT, fT, xT, iT0, iT1, pdT, tT, methodT,Nothing}(sde, f, xs, idx₀, idx₁, pdgrid, t₀, t₁, method, nothing)
end

# Evaluating the intgarals
get_IK_weights!(IK::IntegrationKernel{sdeT}) where sdeT<:AbstractSDE{N,N} where N = quadgk!(IK,IK.temp,IK.xs...)

##############################################
## IK calls
# Scalar problem
function (IK::IntegrationKernel{sdeT})(vals,x₀) where sdeT<:AbstractSDE{1,1}
    basefun_vals!(IK.xs.itp,vals,IK.xs,x₀)
    vals .*= _tp(IK.sde,IK.pdgrid.xs[1][IK.idx₁...],IK.t₁, x₀,IK.t₀, method = IK.method)
    vals
end

function (IK::IntegrationKernel{sdeT})(x₀) where sdeT<:AbstractSDE{1,1}
    basefun_vals(IK.xs.itp,IK.xs,x₀) .* _tp(IK.sde,IK.pdgrid.xs[1][IK.idx₁...],IK.t₁, x₀,IK.t₀, method = IK.method)
end

# 1D Oscillator problem
function (IK::IntegrationKernel{sdeT})(vals,v₀) where sdeT<:SDE_Oscillator1D
    ξ = get_ξ(IK.method,IK.sde,IK.t₁,IK.t₀,IK.pdgrid.xs[1][IK.idx₁[1]],nothing,nothing,v₀) # v₁, x₀ = nothing
    basefun_vals!(IK.pdgrid.xs[1].itp,IK.pdgrid.xs[1].itp.tmp,IK.pdgrid.xs[1], ξ)
    basefun_vals!(IK.pdgrid.xs[2].itp,IK.pdgrid.xs[2].itp.tmp,IK.pdgrid.xs[2],v₀)
    for j in eachindex(IK.pdgrid.xs[2].itp.tmp)
        for i in eachindex(IK.pdgrid.xs[1].itp.tmp)
            vals[i,j] = IK.pdgrid.xs[1].itp.tmp[i] * IK.pdgrid.xs[2].itp.tmp[j] * _tp(IK.sde,IK.sde.par,IK.pdgrid.xs[1][IK.idx₁[1]], IK.pdgrid.xs[2][IK.idx₁[2]], IK.t₁,ξ,v₀,IK.t₀, method = IK.method) 
        end
    end
    
    vals
end

function (IK::IntegrationKernel{sdeT})(v₀) where sdeT<:SDE_Oscillator1D
    basefun_vals(IK.xs.itp,IK.xs,v₀) .* _tp(IK.sde, IK.pdgrid.xs[1][IK.idx₁...], IK.t₁, v₀, IK.t₀, method = IK.method)
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
function _integrate(p::Vp,xs::Tuple) where {Vp<:AbstractArray{Tp,N},xsT<:Axis} where {N,Tp<:Number}
    sp = size(p,N);
    sum(last(xs).wts[i]*_integrate(view(p,(Colon() for _ in 1:N-1)...,i),xs[1:N-1]) for i in 1:sp)
end


