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

struct SDE_VIO{wn,sdeT,wT,idT} <: AbstractSDE{2,2,1} # 1 DoF vibroimpact oscillator
    sde::sdeT
    wall::wT
    ID::idT
    # wT<:Tuple{Wall} = it is assumed, that it is a bottom wall (going down)
end

struct Wall{rT,pT} <: Function
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
"""
    Euler()

Explicit Euler approximation of the drift in a time step: ``x₁ = x₀ + f(x₀, t₀) Δt``.

Pass it as the time stepping `method` of [`PathIntegration`](@ref). See also [`RK2`](@ref), [`RK4`](@ref).
"""
struct Euler <: ExplicitDriftMethod end
"""
    RungeKutta{order}

Explicit Runge–Kutta approximation of the drift in a time step, defined by a Butcher tableau. Construct it with [`RK2`](@ref) or [`RK4`](@ref).

The object also holds the buffers of the stages (resized to the dimension of the SDE): do not share one object between
[`PathIntegration`](@ref)s of different dimensions, or between `PathIntegration`s that are computed at the same time (e.g. in different threads).
"""
struct RungeKutta{order,btT,ksT,tT} <: ExplicitDriftMethod
    BT::btT
    ks::ksT
    temp::tT
end
struct ButcherTableau{aT,bT,cT,_cT}
    a::aT
    b::bT
    c::cT
    _c::_cT
end
struct BTElement{iT,_wT,wT,vT}
    idx::iT
    _weight::_wT
    weight::wT
    val::vT
end

"""
    Maruyama()

Euler–Maruyama approximation of the diffusion in a time step: the noise increment of the last coordinate is ``g(x₀, t₀) ΔW`` with ``ΔW ∼ N(0, Δt)``,
so the transitional PDF is a Gaussian with variance ``g(x₀, t₀)² Δt`` around the result of the drift step ([`Euler`](@ref), [`RK2`](@ref), [`RK4`](@ref)).
It is the (currently only) diffusion approximation, and it is used by default.
"""
struct Maruyama <: ExplicitDiffusionMethod end
struct Milstein <: ExplicitDiffusionMethod end
struct DiscreteTimeStepping{TDrift,TDiff} <: DiscreteTimeSteppingMethod
    drift::TDrift
    diffusion::TDiff
end
abstract type AbstractSDEStep{d,k,m} end
struct NonSmoothSDEStep{d,k,m,sdeT,n,snsT,iddT,idaT,qT} <: AbstractSDEStep{d,k,m}
    sde::sdeT
    sdesteps::snsT
    ID_dyn::iddT
    ID_aux::idaT
    Q_aux::qT
end
struct SDEStep{d,k,m,sdeT,methodT,tracerT,x0T,x1T,tT,tiT,xiT,xi2T} <: AbstractSDEStep{d,k,m}
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
    xi2::xi2T
end

abstract type PreComputeLevel end
struct PreComputeJacobian <: PreComputeLevel end
struct PreComputeLU <: PreComputeLevel end
struct PreComputeNewtonStep <: PreComputeLevel end
abstract type AbstractSymbolicNewtonStepTracer end
struct SymbolicNewtonStepTracer{xIT,xT,detJIT,tempIT,tempT} <: AbstractSymbolicNewtonStepTracer
    xI_0!::xIT
    x_0!::xT
    detJI_inv::detJIT
    tempI::tempIT
    temp::tempT
end
struct VIO_SymbolicNewtonImpactStepTracer{xIT,xT,detJIT,tempIT,tempT,vi_fT,vi_tT,vi_T} <: AbstractSymbolicNewtonStepTracer
    xI_0!::xIT
    x_0!::xT
    detJI_inv::detJIT
    tempI::tempIT
    temp::tempT
    v_toimpact!::vi_fT
    # v_afterimpact!::vs_aT
    vitemp::vi_tT
    v_i::vi_T # [v_beforeimpact, v_afterimpact, v1]
    # v_a::vs_atT
end

# struct StepJacobianLU{JT, JMT}

# end
struct StepJacobian{JT,JMT,tempT}
    J!::JT
    JM::JMT
    temp::tempT
end

# Structs needed for PDF interpolation
abstract type AbstractAxisGrid{T} <: AbstractVector{T} end
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

struct SparseInterpolationBaseVals{Order,vT,iT,lT}
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
struct InterpolatedFunction{T,N,itp_type,axesT,pT,idx_itT,val_itT} <: Function #<:AbstractArray{T,N}
    axes::axesT
    p::pT
    idx_it::idx_itT
    val_it::val_itT
end
"""
    PathIntegration

The response PDF of a stochastic dynamical system together with the step matrices that advance it in time
(step matrix multiplication path integration, see Sykora, Kuske & Yurchenko, Computers and Structures 273 (2022) 106896).
Construct it with `PathIntegration(sde, method, ts, axes...)`.

# Fields
- `pdf`: the response PDF ``p(x, t)`` at `t = PI.t`, an [`InterpolatedFunction`](@ref); `PI(x...)` evaluates it
- `t`: the time
- `ts`: the time points of the step matrices
- `stepMX`: the step matrices, one for each time interval of `ts` (see [`stepMX`](@ref))
- `step_idx`: the index of the step matrix used in the last step
- `marginal_pdfs`: the marginal PDFs (see [`update_mPDFs!`](@ref)), `nothing` if none are requested
- `stepMX_wts`: ``Sᵀw`` for each step matrix ``S`` (``w``: the quadrature weights of the grid), used to normalise the PDF in [`advance!`](@ref)
- `step_dynamics`, `IK`: the time step of the SDE ([`SDEStep`](@ref)) and the integration kernel used to compute the step matrices
- `p_temp`, `kwargs`: a buffer and the keyword arguments of the construction

See also [`advance!`](@ref), [`advance_till_converged!`](@ref), [`steady_state!`](@ref), [`recompute_PI!`](@ref).
"""
mutable struct PathIntegration{dynT,pdT,tsT,stepmxT,Tstp_idx,IKT,ptempT,mpdtT,kwargT,TT,wST}
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
    stepMX_wts::wST # Sᵀw for each step matrix (w: quadrature weights): ∫(S p) = dot(Sᵀw, p)
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
struct IntegrationKernel{kd,sdeT,x1T,diT,fT,pdfT,tT,tempT,kwargT}
    sdestep::sdeT
    x1::x1T
    f::fT # function to integrate over
    discreteintegrator::diT
    t::tT
    pdf::pdfT
    temp::tempT
    kwargs::kwargT
end
struct IK_temp{VT,MT,idxT,valT,kT}
    idx_it::idxT# = Base.Iterators.product(eachindex.(IK.temp.itpVs)...)
    val_it::valT# = Base.Iterators.product(eachindex.(IK.temp.itpVs)...)
    itpVs::VT
    itpM::MT # the row of the step matrix
    kernel::kT # AbstractRowKernel: how the row is integrated
end
# Row kernels: how a row of the step matrix is integrated
abstract type AbstractRowKernel end
# Integrates the full row (length N^d) with the discrete integrator
struct GenericRowKernel <: AbstractRowKernel end
# Sparse interpolations: only the elements touched by the interpolation stencils are accumulated in IK_temp.itpM
struct SparseAccumulator{mT,tT} <: AbstractRowKernel
    mark::mT # mark[j]: element j is touched in the current row
    touched::tT # linear indices of the touched elements
end
# Dense interpolations: the row is the tensor contraction of the basis function values at the quadrature nodes
struct DenseTensorKernel{BT,cT,B1cT,KRT} <: AbstractRowKernel
    Bs::BT # Bs[j][:,q]: basis function values along axis j at the q-th quadrature node
    c::cT # c[q]: quadrature weight × transitional PDF at the q-th quadrature node
    B1c::B1cT # Bs[1] .* transpose(c)
    KR::KRT # Khatri–Rao product of Bs[2:d] (d ≥ 3), nothing otherwise
end
struct Slicer{n,N,idT,slT}
    slicer::slT
end
# struct ImpactInterval{limT,wT}
#     lims::limT # r
#     wallID::wT
#     Q_atwall::BitArray{1}
# end
abstract type AbstractDiscreteIntegratorMethod{dim} end
abstract type AbstractDiscreteIntegratorType{dim} end
"""
    ClenshawCurtisIntegrator(N = 31)

Clenshaw–Curtis quadrature with `N` Chebyshev nodes (the extrema of the Chebyshev polynomial of degree `N - 1`, including both ends of the interval).
Use it as the `discreteintegrator` of [`PathIntegration`](@ref).
"""
struct ClenshawCurtisIntegrator{dim,NT} <: AbstractDiscreteIntegratorMethod{dim}
    N::NT
    function ClenshawCurtisIntegrator(N=31; dim = length(N))
        _N = get_N(N, dim)
        new{dim,typeof(_N)}(_N)
    end
end
"""
    GaussLegendreIntegrator(N = 31)

Gauss–Legendre quadrature with `N` nodes (the ends of the interval are not nodes).
It is the default `discreteintegrator` of [`PathIntegration`](@ref) (with `N = di_N`).
"""
struct GaussLegendreIntegrator{dim,NT} <: AbstractDiscreteIntegratorMethod{dim}
    N::NT
    function GaussLegendreIntegrator(N=31; dim = length(N))
        _N = get_N(N, dim)
        new{dim,typeof(_N)}(_N)
    end
end
"""
    GaussRadauIntegrator(N = 31)

Gauss–Radau quadrature with `N` nodes, including the start of the interval.
Use it as the `discreteintegrator` of [`PathIntegration`](@ref).
"""
struct GaussRadauIntegrator{dim,NT} <: AbstractDiscreteIntegratorMethod{dim}
    N::NT
    function GaussRadauIntegrator(N=31; dim = length(N))
        _N = get_N(N, dim)
        new{dim,typeof(_N)}(_N)
    end
end
"""
    GaussLobattoIntegrator(N = 31)

Gauss–Lobatto quadrature with `N` nodes, including both ends of the interval.
Use it as the `discreteintegrator` of [`PathIntegration`](@ref).
"""
struct GaussLobattoIntegrator{dim,NT} <: AbstractDiscreteIntegratorMethod{dim}
    N::NT
    function GaussLobattoIntegrator(N=31; dim = length(N))
        _N = get_N(N, dim)
        new{dim,typeof(_N)}(_N)
    end
end
"""
    TrapezoidalIntegrator(N = 31)

Composite trapezoidal rule with `N` equidistant nodes.
Use it as the `discreteintegrator` of [`PathIntegration`](@ref).
"""
struct TrapezoidalIntegrator{dim,NT} <: AbstractDiscreteIntegratorMethod{dim}
    N::NT
    function TrapezoidalIntegrator(N=31; dim = length(N))
        _N = get_N(N, dim)
        new{dim,typeof(_N)}(_N)
    end
end
"""
    NewtonCotesIntegrator(N = 31, order = 2)

Composite Newton–Cotes rule with `N` equidistant nodes: `order = 1` is the trapezoidal rule, `2` Simpson's rule, `3` Simpson's 3/8 rule
(with end corrections if `N - 1` is not divisible by `order`). Use it as the `discreteintegrator` of [`PathIntegration`](@ref).
"""
struct NewtonCotesIntegrator{dim,ord,NT} <: AbstractDiscreteIntegratorMethod{dim}
    N::NT
    function NewtonCotesIntegrator(N=31, ord=2; dim = length(N))
        _N = get_N(N, dim)
        new{dim,ord,typeof(_N)}(_N)
    end
end
struct QuadGKIntegrator{iT,rT,kT,qT} <: AbstractDiscreteIntegratorType{1}
    int_limits::iT
    res::rT
    kwargs::kT
    Q_integrate::qT
    res0::rT
end
struct DiscreteIntegrator{dim,xT,wT,resT,tempT,qT} <: AbstractDiscreteIntegratorType{dim}
    x::xT
    w::wT
    res::resT
    temp::tempT
    Q_integrate::qT
    # nodes and weights on the initial interval: rescaling always starts from these
    x_ref::xT
    w_ref::wT
end
struct NonSmoothDiscreteIntegrator{dim,NoDyn,disT} <: AbstractDiscreteIntegratorType{dim}
    discreteintegrators::disT
end

struct DiagonalNormalPDF{uT,sT} <: Function
    μ::uT
    σ²::sT
end

# Step matrix representation types
abstract type StepMatrixRepresentation end
"""
    DenseMX()

Dense step matrix representation, the default for dense interpolations ([`ChebyshevAxis`](@ref), [`TrigonometricAxis`](@ref)) if `d ≤ 2`.
The step matrix ``S`` is stored as `transpose(Sᵀ)` with a dense `Sᵀ` (the columns of `Sᵀ` are the rows of ``S``).
It needs ``8N²`` bytes for ``N`` grid points. Pass it as the `stepMXtype` of [`PathIntegration`](@ref).
"""
struct DenseMX <: StepMatrixRepresentation
end

struct SparseMX{tf,tolT,Ti} <: StepMatrixRepresentation
    Q_threaded::Bool
    tol::tolT # absolute tolerance
    rtol::tolT # tolerance relative to max|S|
end

# Workspaces for computing the rows of the step matrix in blocks of consecutive rows (chunks):
# the rows chunks[c] are computed with IKs[c] (IKs[1] is the original integration kernel), and for a sparse
# step matrix they are stored in buffers[c]
struct RowWorkspace{IKT,bT}
    IKs::Vector{IKT}
    chunks::Vector{UnitRange{Int}}
    buffers::bT
end

# How the rows of the step matrix are computed
abstract type RowComputation end
"""
    SerialRowComputation()

Compute the rows of the step matrix one after the other.
"""
struct SerialRowComputation <: RowComputation end
"""
    ThreadedRowComputation(N_threads = Threads.nthreads(); single_threaded_blas = false)
    ThreadedRowComputation(; N_threads = Threads.nthreads(), single_threaded_blas = false)

Split the rows of the step matrix into `N_threads` blocks of consecutive rows, and compute the blocks in parallel with `Threads.@threads :static`.
Every block has its own copy of the integration kernel (the buffers of the time stepping, the interpolation and the integration).
If it is called inside another `Threads.@threads` loop, the blocks are computed one after the other (with the same result).

With dense interpolations (Chebyshev, trigonometric) every row uses a small BLAS matrix product, and the BLAS threads compete with the Julia threads.
`single_threaded_blas = true` sets `BLAS.set_num_threads(1)` while the blocks are computed in parallel, and restores the number of BLAS threads afterwards
(about 1.5–2× faster with 16 threads). The number of BLAS threads is a global setting: do not use it when BLAS is used concurrently
(e.g. several `PathIntegration`s are computed in parallel with `Threads.@spawn`), as the setting could be restored incorrectly.
"""
struct ThreadedRowComputation <: RowComputation
    N_threads::Int
    single_threaded_blas::Bool
    function ThreadedRowComputation(N::Integer = Threads.nthreads(); N_threads::Integer = N, single_threaded_blas::Bool = false)
        N_threads ≥ 1 || throw(ArgumentError("N_threads = $N_threads, it has to be at least 1"))
        new(N_threads, single_threaded_blas)
    end
end
"""
    BatchRowComputation(N_threads = Threads.nthreads(); single_threaded_blas = false)
    BatchRowComputation(; N_threads = Threads.nthreads(), single_threaded_blas = false)

Split the rows of the step matrix into `N_threads` blocks of consecutive rows, and compute the blocks in parallel with `Polyester.@batch`.
Every block has its own copy of the integration kernel (the buffers of the time stepping, the interpolation and the integration).
Polyester only uses the threads that are free, e.g. inside another threaded loop the blocks are computed one after the other (with the same result).

`single_threaded_blas`: see [`ThreadedRowComputation`](@ref).
"""
struct BatchRowComputation <: RowComputation
    N_threads::Int
    single_threaded_blas::Bool
    function BatchRowComputation(N::Integer = Threads.nthreads(); N_threads::Integer = N, single_threaded_blas::Bool = false)
        N_threads ≥ 1 || throw(ArgumentError("N_threads = $N_threads, it has to be at least 1"))
        new(N_threads, single_threaded_blas)
    end
end