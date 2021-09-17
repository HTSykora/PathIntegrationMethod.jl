struct Scalar_Or_Function{fT} <: Function
    f::fT
end
TupleVectorUnion = Union{Tuple,AbstractVector}

# Types to describe the SDE dynamics
abstract type AbstractSDE{d,k,m} end
# f: Rᵈ×[0,T] ↦ Rᵈ, g: Rᵈ × [0,T] ↦ Rᵈˣᵐ, gᵢ,ⱼ = 0 for i = 1,...,m; j = k ... d
struct SDE{d,k,m,fT,gT,pT} <: AbstractSDE{d,k,m}
    f::fT
    g::gT
    par::pT
end

struct SDE_Oscillator1D{fT, gT, parT} <: AbstractSDE{2,2,1}
    f::fT
    g::gT
    par::parT
end

struct SDE_VI_Oscillator1D{wT, oscT} <: AbstractSDE{2,2,1}
    osc1D::oscT
    wall::wT
    # wT<:Tuple{Wall} = it is assumed, that it is a bottom wall (going down)
end
struct Wall{rT,pT}
    r::rT
    pos::pT
    # impact_v_sign::dT
end
abstract type AbstractSDEComponent end
struct DriftTerm{d,fT} <: AbstractSDEComponent
    f::fT
end
struct DiffusionTerm{d,k,kd,m,gT} <: AbstractSDEComponent
    g::gT
end

abstract type DiscreteTimeSteppingMethod end
abstract type DiscreteDriftMethod <: DiscreteTimeSteppingMethod end
abstract type DiscreteDiffusionMethod <: DiscreteTimeSteppingMethod end
struct Euler <: DiscreteDriftMethod end
struct RK4 <: DiscreteDriftMethod end
struct Maruyama <: DiscreteDiffusionMethod end
struct Milstein <: DiscreteDiffusionMethod end
struct DiscreteTimeStepping{TDrift,TDiff} <:DiscreteTimeSteppingMethod
    drift::TDrift
    diffusion::TDiff
end

struct SDEStep{d, k, m, sdeT, methodT,tracerT,x0T,x1T,tT}
    sde::sdeT
    method::methodT
    x0::x0T
    x1::x1T
    t0::tT
    t1::tT
    
    steptracer::tracerT # Required for discrete time step backtracing
end
abstract type PreComputeLevel end
struct PreComputeJacobian <: PreComputeLevel end
struct PreComputeLU <: PreComputeLevel end
struct PreComputeNewtonStep <: PreComputeLevel end
struct SymbolicNewtonStep{xIT, detJiT,tempT}
    xI_0!::xIT
    # xII_1!::xIIT
    detJ_inv::detJiT
    temp::tempT
end
# struct StepJacobianLU{JT, JMT}
    
# end
struct StepJacobian{JT, JMT,tempT}
    J!::JT
    JM::JMT
    temp::tempT
end

# Structs needed for PDF interpolation
abstract type AbstractGridAxis{T} <:AbstractVector{T} end 
struct GridAxis{itpT,wT,xT,xeT,tmpT} <: AbstractGridAxis{xeT}
    itp::itpT
    xs::xT
    wts::wT
    temp::tmpT
end

abstract type AbstractInterpolationType end
struct ChebyshevInterpolation{N} <: AbstractInterpolationType end
struct FourierInterpolation{N} <: AbstractInterpolationType end
struct LinearInterpolation{ΔT} <: AbstractInterpolationType
    Δ::ΔT
end
struct TrapezoidalWeights{T,ΔT} <: AbstractVector{T}
    l::Int64
    Δ::ΔT
end

struct ProbabilityDensityFunction{T,N,axesT,pT,idx_itT} <: AbstractArray{T,N}
    axes::axesT
    p::pT
    idx_it::idx_itT
end
struct PathIntegrationProblem{dynT, pdT, tsT, tpdMX_type, Tstp_idx}
    step_dynamics::dynT # SDEStep
    pdgrid::pdT
    ts::tsT
    step_MX::tpdMX_type
    step_idx::Tstp_idx
end

# Utility types
struct IntegrationKernel{sdeT,xT,fT,pdfT,x1idxT,tidxT,tempT,iiT}
    sdestep::sdeT
    f::fT # function to integrate over
    xs::xT # integration axes
    x1_idx::x1idxT
    t_idx::tidxT
    pdf::pdfT
    temp::tempT
end
struct IK_temp{VT,MT}
    itpVs::VT
    itpM::MT
    # impactinterval::iiT
end
# struct ImpactInterval{limT,wT}
#     lims::limT # r
#     wallID::wT
#     Q_atwall::BitArray{1}
# end



