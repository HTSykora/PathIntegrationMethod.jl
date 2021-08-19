@inline function cleanup_quadgk_keywords(;σ_init = nothing, kwargs...)
    (;kwargs...)
end

# Integration kernel for custom kernel functions
function IntegrationKernel(f::fT,xs, tempsize; kwargs...) where fT<: Function
    temp = zeros(Float64,tempsize...) # get eltype properly
    IntegrationKernel(nothing,f,xs,nothing,nothing,nothing,nothing,nothing,nothing,temp, nothing)
end

# Integration kernel for fast computation of the transition matrix
function IntegrationKernel(sde::sdeT,f::fT,xs::xT,idx₀::iT0,idx₁::iT1,pdgrid::pdT,t₀::tT,t₁::tT,method::methodT; kwargs...) where {sdeT, fT, xT, iT0, iT1, pdT, tT, methodT}
    IntegrationKernel{sdeT, fT, xT, iT0, iT1, pdT, tT, methodT,Nothing,Nothing}(sde, f, xs, idx₀, idx₁, pdgrid, t₀, t₁, method, nothing,nothing)
end

@inline get_t0(IK::IntegrationKernel) = IK.t₀
@inline get_t1(IK::IntegrationKernel) = IK.t₁

@inline get_t0(IK::IntegrationKernel{sdeT,iT0, iT1,xT,fT,pdT,tT}) where {sdeT,iT0, iT1,xT,fT,pdT, tT<: Vector{teT}} where teT<:Number = IK.t₀[1]
@inline get_t1(IK::IntegrationKernel{sdeT,iT0, iT1,xT,fT,pdT,tT}) where {sdeT,iT0, iT1,xT,fT,pdT, tT<: Vector{teT}} where teT<:Number = IK.t₁[1]

# Evaluating the integrals
function get_IK_weights!(IK::IntegrationKernel{sdeT}; integ_limits = (IK.xs[1],IK.xs[end]), kwargs...) where sdeT<:AbstractSDE{N,N} where N 
    # integ_limits = IK.xs -> # old version
    quadgk!(IK,IK.temp,integ_limits...; cleanup_quadgk_keywords(;kwargs...)...)
end
##############################################
## IK calls
# Generic functions
function (IK::IntegrationKernel{Nothing,Nothing,Nothing,xT,fT})(vals,xs...) where {xT,fT<: Function}
    for (i,x) in enumerate(xs)
        basefun_vals_safe!(IK.xs[i].itp,IK.xs[i].itp.tmp,IK.xs[i],x)
    end
    idx_it = Base.Iterators.product(eachindex.(IK.xs)...)
    fx = IK.f(xs...);

    for (i,idx) in enumerate(idx_it)
        vals[idx...] = fx*prod(IK.xs[j].itp.tmp[id] for (j,id) in enumerate(idx))
    end
    vals
end

# Scalar problem
function (IK::IntegrationKernel{sdeT})(vals,x₀) where sdeT<:AbstractSDE{1,1}
    basefun_vals!(IK.xs.itp,vals,IK.xs,x₀)
    fx = _tp(IK.sde,IK.pdgrid.xs[1][IK.idx₁...],get_t1(IK), x₀,get_t0(IK), method = IK.method);
    vals .*= fx
    vals
end

# 1D Oscillator problem + vibroimpact
function (IK::IntegrationKernel{sdeT})(vals,v₀) where sdeT<:Union{SDE_Oscillator1D}
    ξ, extra_args... = get_ξ(IK.method,IK.sde,get_t1(IK),get_t0(IK),IK.pdgrid.xs[1][IK.idx₁[1]],nothing,nothing,v₀) # v₁, x₀ = nothing
    # for impact system:
    # extra_args  == Q_impact, Δt1, Δt2, r

    basefun_vals_safe!(IK.pdgrid.xs[1].itp,IK.pdgrid.xs[1].itp.tmp,IK.pdgrid.xs[1], ξ)
    basefun_vals!(IK.pdgrid.xs[2].itp,IK.pdgrid.xs[2].itp.tmp,IK.pdgrid.xs[2],v₀)

    fx = _tp(IK.sde,_par(IK.sde),IK.pdgrid.xs[1][IK.idx₁[1]], IK.pdgrid.xs[2][IK.idx₁[2]], get_t1(IK),ξ,v₀,get_t0(IK), extra_args...; method = IK.method) 
    # for impact system:
    # fx = _tp(IK.sde,_par(IK.sde),IK.pdgrid.xs[1][IK.idx₁[1]], IK.pdgrid.xs[2][IK.idx₁[2]], get_t1(IK),ξ,v₀,get_t0(IK), Q_impact, Δt1, Δt2, r , method = IK.method) 

    for j in eachindex(IK.pdgrid.xs[2].itp.tmp)
        for i in eachindex(IK.pdgrid.xs[1].itp.tmp)
            vals[i,j] = IK.pdgrid.xs[1].itp.tmp[i] * IK.pdgrid.xs[2].itp.tmp[j] * fx
        end
    end
    
    vals
end

# 1D VibroImpact Oscillator problem
function (IK::IntegrationKernel{sdeT})(vals,v₀) where sdeT<:Union{SDE_VI_Oscillator1D}
    ξ = get_ξ(IK.method,IK.sde.osc1D, get_t1(IK),get_t0(IK),IK.pdgrid.xs[1][IK.idx₁[1]],nothing,nothing,v₀) # v₁, x₀ = nothing
    # for impact system:
    # extra_args  == Q_impact, Δt1, Δt2, r

    basefun_vals_safe!(IK.pdgrid.xs[1].itp,IK.pdgrid.xs[1].itp.tmp,IK.pdgrid.xs[1], ξ)
    basefun_vals!(IK.pdgrid.xs[2].itp,IK.pdgrid.xs[2].itp.tmp,IK.pdgrid.xs[2],v₀)

    fx = _tp(IK.sde,_par(IK.sde),IK.pdgrid.xs[1][IK.idx₁[1]], IK.pdgrid.xs[2][IK.idx₁[2]], get_t1(IK),ξ,v₀,get_t0(IK); method = IK.method) 
    # for impact system:
    # fx = _tp(IK.sde,_par(IK.sde),IK.pdgrid.xs[1][IK.idx₁[1]], IK.pdgrid.xs[2][IK.idx₁[2]], get_t1(IK),ξ,v₀,get_t0(IK), Q_impact, Δt1, Δt2, r , method = IK.method) 
    for j in eachindex(IK.pdgrid.xs[2].itp.tmp)
        for i in eachindex(IK.pdgrid.xs[1].itp.tmp)
            vals[i,j] = IK.pdgrid.xs[1].itp.tmp[i] * IK.pdgrid.xs[2].itp.tmp[j] * fx
        end
    end
    
    if IK.wallID[1] != 0
        ξ, Δt₁, Δt₂, r = get_ξ_impact(IK.method,IK.sde.osc1D, get_t1(IK),get_t0(IK),IK.pdgrid.xs[1][IK.idx₁[1]],nothing,nothing,v₀,IK.wallID[1]) # v₁, x₀ = nothing
        basefun_vals_safe!(IK.pdgrid.xs[1].itp,IK.pdgrid.xs[1].itp.tmp,IK.pdgrid.xs[1], ξ)
        fx = _tp(IK.sde,_par(IK.sde),IK.pdgrid.xs[1][IK.idx₁[1]], IK.pdgrid.xs[2][IK.idx₁[2]], get_t1(IK),ξ,v₀,get_t0(IK), Δt₁, Δt₂, r; method = IK.method) / r(v₀)

        for j in eachindex(IK.pdgrid.xs[2].itp.tmp)
            for i in eachindex(IK.pdgrid.xs[1].itp.tmp)
                vals[i,j] = vals[i,j] + IK.pdgrid.xs[1].itp.tmp[i] * IK.pdgrid.xs[2].itp.tmp[j] * fx
            end
        end
    end
    vals
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

# Integrate over a full PDGrid
function _integrate(p::PDGrid)
    _integrate(p.p,p.xs)
end
