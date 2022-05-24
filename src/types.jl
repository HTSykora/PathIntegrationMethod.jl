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

# struct SDE_Oscillator1D{fT, gT, parT} <: AbstractSDE{2,2,1}
#     f::fT
#     g::gT
#     par::parT
# end

struct SDE_VIO{sdeT, wT} <: AbstractSDE{2,2,1} # 1 DoF vibroimpact oscillator
    sde::sdeT
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
abstract type ExplicitDriftMethod <: DiscreteTimeSteppingMethod end
abstract type ExplicitDiffusionMethod <: DiscreteTimeSteppingMethod end
struct Euler <: ExplicitDriftMethod end
struct RungeKutta{order,btT,ksT,tT} <: ExplicitDriftMethod
    BT::btT
    ks::ksT
    temp::tT
end
struct ButcherTableau{aT, bT, cT,_cT}
    a::aT
    b::bT
    c::cT
    _c::_cT
end
struct BTElement{iT,_wT, wT,vT}
    idx::iT
    _weight::_wT
    weight::wT
    val::vT
end

struct Maruyama <: ExplicitDiffusionMethod end
struct Milstein <: ExplicitDiffusionMethod end
struct DiscreteTimeStepping{TDrift,TDiff} <:DiscreteTimeSteppingMethod
    drift::TDrift
    diffusion::TDiff
end
struct NonSmoothSDEStep{d,k,m,snsT}
    sdesteps::snsT
end
struct SDEStep{d, k, m, sdeT, methodT,tracerT,x0T,x1T,tT, tiT, xiT}
    sde::sdeT
    method::methodT
    x0::x0T
    x1::x1T
    t0::tT
    t1::tT
    
    steptracer::tracerT # Required for discrete time step backtracing
    # intermediate step utilities
    ti::tiT
    xi::xiT
end
abstract type PreComputeLevel end
struct PreComputeJacobian <: PreComputeLevel end
struct PreComputeLU <: PreComputeLevel end
struct PreComputeNewtonStep <: PreComputeLevel end
struct SymbolicNewtonStepTracer{xIT, xT, detJIT,tempIT,tempT}
    xI_0!::xIT
    x_0!::xT
    # xII_1!::xIIT
    detJI_inv::detJIT
    tempI::tempIT
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
abstract type AbstractAxisGrid{T} <:AbstractVector{T} end 
struct AxisGrid{itpT,wT,xT,xeT,tmpT} <: AbstractAxisGrid{xeT}
    itp::itpT
    xs::xT
    wts::wT
    temp::tmpT
end

abstract type AbstractInterpolationType end
abstract type DenseInterpolationType <: AbstractInterpolationType end
abstract type SparseInterpolationType <: AbstractInterpolationType end
struct ChebyshevInterpolation{NT,_1T} <: DenseInterpolationType 
    N::NT
    _1::_1T
end
struct TrigonometricInterpolation{_isodd,NT,cT} <: DenseInterpolationType
    N::NT
    c::cT
end
struct LinearInterpolation{ΔT} <: SparseInterpolationType
    Δ::ΔT
end
# struct TrapezoidalWeights{ΔT} <: AbstractVector{ΔT}
#     l::Int64
#     Δ::ΔT
# end
struct NewtonCotesWeights{N,ΔT,remainder,lT} <: AbstractVector{ΔT}
    # remainder: mod(l-1,N) can be handled by only a pure Newton Cotes formula (e.g. compatible length with the integration scheme)
    l::lT
    Δ::ΔT
end

struct SparseInterpolationBaseVals{Order,vT, iT, lT}
    val::vT
    idxs::iT
    l::lT
end
struct CubicInterpolation{ΔT} <: SparseInterpolationType
    Δ::ΔT
end
struct QuinticInterpolation{ΔT} <: SparseInterpolationType
    Δ::ΔT
end
struct InterpolatedFunction{T,N,itp_type,axesT,pT,idx_itT,val_itT} <:Function #<:AbstractArray{T,N}
    axes::axesT
    p::pT
    idx_it::idx_itT
    val_it::val_itT
end
mutable struct PathIntegration{dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT,kwargT,TT}
    step_dynamics::dynT # SDEStep
    pdf::pdT
    p_temp::ptempT
    ts::tsT
    stepMX::stepmxT
    step_idx::Tstp_idx
    IK::IKT
    marginal_pdfs::mpdtT
    kwargs::kwargT
    t::TT
end

struct MarginalPDF{pT,idT,wT,tT,p0T,dT}
    pdf::pT
    ID::idT
    wMX::wT
    temp::tT
    p0::p0T
    dims::dT
end
# Utility types
struct IntegrationKernel{kd,sdeT,x1T,diT,fT,pdfT, tT,tempT,kwargT}
    sdestep::sdeT
    x1::x1T
    f::fT # function to integrate over
    discreteintegrator::diT
    t::tT
    pdf::pdfT
    temp::tempT
    kwargs::kwargT
end
struct IK_temp{VT,MT,idxT,valT}
    idx_it::idxT# = Base.Iterators.product(eachindex.(IK.temp.itpVs)...)
    val_it::valT# = Base.Iterators.product(eachindex.(IK.temp.itpVs)...)
    itpVs::VT
    itpM::MT
    # impactinterval::iiT
end
struct Slicer{n,N,idT,slT}
    slicer::slT
end
# struct ImpactInterval{limT,wT}
#     lims::limT # r
#     wallID::wT
#     Q_atwall::BitArray{1}
# end
abstract type AbstractDiscreteIntegratorMethod end
abstract type AbstractDiscreteIntegratorType end
struct ClenshawCurtisIntegrator <:AbstractDiscreteIntegratorMethod end
struct GaussLegendreIntegrator <:AbstractDiscreteIntegratorMethod end
struct GaussRadauIntegrator <:AbstractDiscreteIntegratorMethod end
struct GaussLobattoIntegrator <:AbstractDiscreteIntegratorMethod end
struct TrapezoidalIntegrator <: AbstractDiscreteIntegratorMethod end
struct NewtonCotesIntegrator{N} <: AbstractDiscreteIntegratorMethod end
struct QuadGKIntegrator{iT,rT,kT} <:AbstractDiscreteIntegratorType
    int_limits::iT
    res::rT
    kwargs::kT
end
struct DiscreteIntegrator{IType,xT,wT,resT,tempT} <:AbstractDiscreteIntegratorType
    x::xT
    w::wT
    res::resT
    temp::tempT
end
struct DiagonalNormalPDF{uT, sT} <: Function
    μ::uT
    σ²::sT
end

# Step matrix representation types
abstract type StepMatrixRepresentation end
struct DenseMX <:StepMatrixRepresentation
end

struct SparseMX{tf,tolT} <:StepMatrixRepresentation
    Q_threaded::Bool
    tol::tolT
end