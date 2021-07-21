struct PathIntegrationProblem{N,k,sdeT,pdT,tpdMX_type,dtT,methodT}
    sde::sdeT
    pdgrid::pdT
    tpdMX::tpdMX_type
    Δt::dtT
    method::methodT
end

struct PDGrid{N,k,T,xT,pT,ξT,gridT,iT} <: AbstractArray{T,N}
    xs::xT
    p::pT
    p_temp::pT
    ξ_temp::ξT
    grid::gridT
    i_temp::iT
end

abstract type AbstractSDE{N,k} end
struct DummySDE{N,k} <: AbstractSDE{N,k} end

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
struct Axis{itpT,xeT,ewT,xT,wT} <: AbstractAxis{xeT}
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

struct IntegrationKernel{sdeT,iT0,iT1,xT,fT,pdT,tT,methodT,tempT}
    sde::sdeT
    f::fT
    xs::xT
    idx₀::iT0
    idx₁::iT1
    pdgrid::pdT
    t₀::tT
    t₁::tT
    method::methodT
    temp::tempT
end

abstract type SDEMethod end
struct EulerMaruyama <: SDEMethod end
struct Milstein <: SDEMethod end
struct RKMaruyama <: SDEMethod end
