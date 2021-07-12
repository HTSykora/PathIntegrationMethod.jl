abstract type AbstractSDE{N,k} end
struct SDE{N,k,fT,gT,pT} <: AbstractSDE{N,k}
    f::fT
    g::gT
    par::pT
end
struct SDE_Oscillator1D{fT, gT, parT} <: AbstractSDE{2,2}
    f::fT
    g::gT
    par::parT
end

abstract type AbstractSDEComponent{N,k,T} end
struct DriftTerm{N,k,fT} <: AbstractSDEComponent{N,k,fT}
    f::fT
end
struct DiffusionTerm{N,k,gT} <: AbstractSDEComponent{N,k,gT}
    g::gT
end

abstract type AbstractAxis{T} <:AbstractVector{T} end 
struct Axis{itpT,xeT,ewT,xT,wT} <: AbstractAxis{T}
    itp::itpT
    x::xT
    wts::wT
end

abstract type AbstractInterpolationType{TF} end # Q_equidistant = TF
struct ChebyshevInterpolation{N,Ttmp} <: AbstractInterpolationType{false}
    tmp::Ttmp
end
struct LinearInterpolation{TF,Ttmp,tΔ} <: AbstractInterpolationType{TF}
    tmp::Ttmp
    Δ::tΔ
end
struct TrapezoidalWeights{T,ΔT} <: AbstractVector{T}
    l::Int64
    Δ::ΔT
end

struct PDGrid{N,k,T,xT,pT,ξT,gridT,iT} <: AbstractArray{T,N}
    xs::xT
    p::pT
    p_temp::pT
    ξ_temp::ξT
    grid::gridT
    i_temp::iT
end

struct IntegrationKernel{iT,xT,fT,ttT,tempT}
    idx₀::iT
    idx₁::iT
    xs::xT
    f::fT
    TT::ttT
    temp::tempT
end

abstract type SDEMethod end
struct EulerMaruyama <: SDEMethod end
struct Milstein <: SDEMethod end
struct RKMaruyama <: SDEMethod end