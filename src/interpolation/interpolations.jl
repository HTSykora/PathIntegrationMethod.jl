abstract type AbstractInterpolationType{TF} end

struct ChebyshevInterpolation{N,Ttmp,wT} <: AbstractInterpolationType{false}
    tmp::Ttmp
    wts::Vector{wT}
end
ChebyshevInterpolation(N) = ChebyshevInterpolation{N-1,Nothing,Float64}(nothing,clenshawcurtisweights(N))
function ChebyshevInterpolation(N,a::T,b::T) where T<:Number
    wts = clenshawcurtisweights(a,b,N)
    ChebyshevInterpolation{N-1,Nothing,eltype(wts)}(nothing,wts)
end

struct TrapezoidalWeights{T} <: AbstractVector{T}
    l::Int64
end
Base.getindex(t::TrapezoidalWeights{T},idx::Integer) where T = (idx == 1 || idx == t.l) ? T(0.5) : T(1.)
Base.size(t::TrapezoidalWeights) = (t.l,)

struct LinearInterpolation{TF,Ttmp,tΔ,wT} <: AbstractInterpolationType{TF}
    tmp::Ttmp
    Δ::tΔ
    wts::TrapezoidalWeights{wT}
end

LinearInterpolation(n; Q_equidistant = false) = LinearInterpolation{Q_equidistant,Nothing,Nothing,Float64}(nothing,nothing,TrapezoidalWeights{Float64}(n))
function LinearInterpolation(n, Δx; Q_equidistant = false) 
    wT = Float64
    _Δx = Q_equidistant ? Δx : nothing
    LinearInterpolation{Q_equidistant,Nothing,typeof(_Δx),wT}(nothing,_Δx,TrapezoidalWeights{wT}(n))
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
            return true, i, zero(Tp)
        end
    else
        return true, 0, zero(Tp)
    end
end


struct IntegrationKernel{iT,xT,fT,ttT,tempT}
    idx₀::iT
    idx₁::iT
    xs::xT
    f::fT
    TT::ttT
    temp::tempT
end

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