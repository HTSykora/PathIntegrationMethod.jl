"""
    PathIntegration(sde, method, ts, axes...; kwargs...)
    PathIntegration(sdestep, ts, axes...; kwargs...)

Set up the computation of the response probability density function (PDF) of the stochastic dynamical system `sde`:
initialise the PDF on the grid spanned by `axes`, and compute the step matrices that advance it in time (see [`advance!`](@ref)).
Returns the [`PathIntegration`](@ref) object `PI`: `PI.pdf` is the response PDF at time `PI.t` (initially 0), and `PI(x...)` evaluates it.

# Arguments
- `sde`: `d`-dimensional stochastic dynamical system, see [`SDE`](@ref) and [`SDE_VIO`](@ref).
- `method`: time stepping used to approximate the transitional PDF: [`Euler`](@ref)`()`, [`RK2`](@ref)`()` or [`RK4`](@ref)`()` for the drift
  (with the [`Maruyama`](@ref) approximation of the diffusion). A method whose order equals the dimension is recommended: `Euler()` for `d = 1`, `RK2()` for `d = 2`, `RK4()` for `d ≥ 3`.
- `ts`: a time step `Δt` (a single step matrix for the time interval `[0, Δt]`), or time points `[t₀, t₁, …, tₙ]` (a step matrix for each interval `[tⱼ₋₁, tⱼ]`,
  which [`advance!`](@ref) uses cyclically, e.g. the time points of one period for time-periodic systems).
- `axes`: `d` [`AxisGrid`](@ref)s (e.g. [`QuinticAxis`](@ref), [`ChebyshevAxis`](@ref)) spanning the region where the PDF is computed.
- `sdestep`: an [`SDEStep`](@ref) (the time stepping of an SDE) instead of `sde` and `method`.

# Keyword Arguments
- `discreteintegrator = defaultdiscreteintegrator(sde, di_N = 31)`: Discrete integrator to evaluate the Chapman-Kolmogorov equation
    - The default discrete integrator algorithms are
        - `defaultdiscreteintegrator(sde::AbstractSDE{d,k,m}, di_N = 31) = GaussLegendreIntegrator(di_N)`
        - `defaultdiscreteintegrator(sde::SDE_VIO, di_N = 31) = Tuple(GaussLegendreIntegrator(di_N) for _ in 1:2)`
    - `di_N = 31`: Resolution of the discrete integrator. Can be a `Integer` or `NTuple{d-k+1,<:Integer}` that defines the discrete integration resolution in each `d-k+1` integration direction.
- `smart_integration = true`: Only integrate where the transitional PDF has nonzero elements. It is approximated with the step function. Use `false` if (time step * diffusion) results in a wide TPDF. Usually `true` is the better choice.
- `int_limit_thickness_multiplier = 6`: The "thickness" scaling of the TPDF during smart integration.
- `initialise_pdf = true`: Initialise the response probability density function (RPDF). If false, then the RPDF is initialised as p(x) ≡ 0.
- `f_init = nothing`: Initial RPDF as a function. If `f_init` is a `Nothing` and `initialise_pdf = true` then a diagonal Gaussian distribution is used as initial distribution
- `μ_init = nothing`: Mean of the initial Gaussian distribution used in case `f_init = nothing`. 
    - `Nothing`: Uses the middle of the range defined by the axis
    - `Number`: Uses the single value `μ_init` for each axis direction
    - `Union{NTuple{d,<:Number},AbstractVector{<:Number}}`: Individual means for each axis.
- `σ_init = nothing`: Standard deviation of the initial diagonal Gaussian distribution used in case `f_init = nothing`. 
    - `Nothing`: Uses the 1/12th of the range width in each direction defined by the `axes`.
    - `Number`: Uses the single value `σ_init` for each axis direction
    - `Union{NTuple{d,<:Number},AbstractVector{<:Number}}`: Individual standard deviations for each axis.
- `pre_compute = true`: Compute the `stepMX`. This should be left unchanged if the RPDF computation is the goal. (With `false`, only the PDF is initialised.)
- `stepMXtype = nothing`: Step matrix representation type. The default representation depends on `d` and the interpolation used (see Interpolation): `d ≤ 2` with sparse interpolation the `stepMX` is a multithreaded sparse matrix, for dense interpolations it is dense, and for `d>2` the default is a multithreaded sparse
    Possible options:
    - `SparseMX(; threaded = true, sparse_tol = 1e-6, sparse_rtol = 0.0, index_type = Int)`
    - `DenseMX()`: dense matrix, stored as `transpose(Sᵀ)`
- `multithreaded_sparse = true`: is the sparse `stepMX` multithreaded if no `stepMXtype` is specified
- `sparse_tol = 1e-6`: absolute tolerance for the elements considered as zero values in the sparse stepMX if no `stepMXtype` is specified
- `sparse_rtol = 0.0`: tolerance relative to the largest element: elements with `|Sᵢⱼ| ≤ max(sparse_tol, sparse_rtol * max|Sᵢⱼ|)` are considered as zero values in the sparse stepMX if no `stepMXtype` is specified. Eq. (37) of Sykora et al. (2022) corresponds to `sparse_tol = 0, sparse_rtol = 1e-8`.
- `index_type = Int`: index type of the sparse stepMX if no `stepMXtype` is specified. `Int32` needs less memory (and memory bandwidth in `advance!`), but it limits the number of nonzero elements to 2³¹-1.
- `rowcomputation = Threads.nthreads() > 1 ? ThreadedRowComputation() : SerialRowComputation()`: how the rows of the stepMX are computed (also used by `recompute_stepMX!`)
    - `SerialRowComputation()`: one row after the other
    - `ThreadedRowComputation(N_threads = Threads.nthreads(); single_threaded_blas = false)`: `N_threads` blocks of rows in parallel using `Threads.@threads`
    - `BatchRowComputation(N_threads = Threads.nthreads(); single_threaded_blas = false)`: `N_threads` blocks of rows in parallel using `Polyester.@batch`
    - `single_threaded_blas = true`: use a single BLAS thread during the parallel computation (faster with dense interpolations, see `ThreadedRowComputation`)
    The result does not depend on the row computation or on the number of threads. (`multithreaded_sparse` controls the multithreading of `advance!`.)
- `backtracing = NewtonBacktracing()`: how the start points of the time steps are traced back from the grid points
  (with an `sdestep`, set it in [`SDEStep`](@ref)`(sde, method, ts; backtracing)`)
    - [`NewtonBacktracing`](@ref)`()`: Newton iteration, compiled from symbolic derivatives
    - [`ExplicitBacktracing`](@ref)`()`: one explicit drift step backward in time, the same transitional PDF (with `RK2()`, `RK4()`), 2–3.5× faster
    - [`StrangSplitting`](@ref)`()`: Strang splitting of the drift and the diffusion with explicit backward drift steps, second order in `Δt` (with `RK2()`, `RK4()` and additive noise)
- `mPDF_IDs = nothing`: marginal PDFs (mPDFs) of the coordinates specified by `mPDF_IDs` (computed by [`update_mPDFs!`](@ref))
    - `Nothing`: no mPDF is initialised
    - `Integer`, e.g. `2`: 1-dimensional mPDF of the coordinate `mPDF_IDs`
    - `NTuple{n,<:Integer}`, e.g. `(1, 2)`: `n`-dimensional mPDF of the coordinates in `mPDF_IDs`
    - a `Vector` (or `Tuple`) of these, e.g. `[1, 2]`: several mPDFs
- `allow_extrapolation::Bool = false`, `zero_extrapolation::Bool = true`: extrapolation of the PDF outside of the region spanned by `axes`
  in the step matrix computation (see [`InterpolatedFunction`](@ref)); by default the PDF is zero outside of the region.

# Example
```julia
f1(x, p, t) = x[2]
f2(x, p, t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x, p, t) = p[3]
sde = SDE((f1, f2), g2, [0.5, 0.25, 1.0])
PI = PathIntegration(sde, RK4(), 0.02, QuinticAxis(-4., 4., 41), QuinticAxis(-4., 4., 41))
advance!(PI)        # one time step
steady_state!(PI)   # the stationary PDF
PI(0.5, 0.0)        # the PDF at x = 0.5, v = 0
```
"""
function PathIntegration(sdestep::AbstractSDEStep{d,k,m}, _ts, axes::Vararg{Any,d}; 
    di_N = 31, discreteintegrator = defaultdiscreteintegrator(sdestep.sde, di_N),
    initialise_pdf = true, f_init = nothing, pre_compute = true, stepMXtype = nothing, sparse_tol = 1e-6, sparse_rtol = 0.0,
    mPDF_IDs = nothing, extract_IK = Val{false}(), rowcomputation = default_rowcomputation(), generic_row_kernel = false, kwargs...) where {d,k,m}
    # (the back-tracing is set by the SDE step)
    if haskey(kwargs, :backtracing) && kwargs[:backtracing] != backtracing(sdestep)
        throw(ArgumentError("backtracing = $(kwargs[:backtracing]), but the SDE step uses $(backtracing(sdestep)): set it in SDEStep(sde, method, ts; backtracing)"))
    end
    if stepMXtype isa StepMatrixRepresentation
        _stepMXtype = stepMXtype
    else
        _stepMXtype = get_stepMXtype(sdestep.sde, get_val_itp_type(axes); sparse_tol = sparse_tol, sparse_rtol = sparse_rtol, kwargs...)
    end

    if initialise_pdf
        if f_init isa Nothing
            _f = init_DiagonalNormalPDF(axes...; kwargs...)
        else
            _f = f_init
        end
    else
        _f = nothing
    end
    # for j in eachindex(sdestep.steptracer.temp)
    #     sdestep.steptracer.temp[j] = zero(sdestep.steptracer.temp[j])
    # end
    # for j in eachindex(sdestep.steptracer.tempI)
    #     sdestep.steptracer.tempI[j] = zero(sdestep.steptracer.tempI[j])
    # end

    pdf = InterpolatedFunction(axes...; f = _f, kwargs...)
    ts = get_ts(_ts);
    step_idx = 0
    t = 0.
    if pre_compute
        IK = integration_kernel(sdestep, discreteintegrator, pdf, ts, generic_row_kernel,
            (;sparse_tol = get_tol(_stepMXtype), sparse_rtol = get_rtol(_stepMXtype), rowcomputation = rowcomputation, kwargs...); kwargs...)
        
        if extract_IK isa Val{true}
            return IK
        end
        stepMX = compute_stepMX(IK; stepMXtype = _stepMXtype, IK.kwargs...)
    else
        stepMX = nothing
        IK = nothing
    end
    if mPDF_IDs isa Nothing
        mpdf = nothing
    else
        mpdf = initialise_mPDF(pdf,mPDF_IDs)
    end

    p_temp = similar(pdf.p);
    PathIntegration(sdestep, pdf, p_temp,ts, stepMX, step_idx, IK, mpdf, kwargs,t, stepMX_weights(stepMX, pdf))
end
_val(vals) = vals
get_ts(_ts::AbstractVector{tsT}) where tsT<:Number = collect(_ts)
get_ts(_ts::tsT) where tsT<:Number = [zero(_ts), _ts]
"""
    stepMX(PI, i = 1)

The step matrix ``S`` of the `i`-th time interval of `PI.ts`. It maps the PDF values at the grid nodes from the start to the end of the interval:
`vec(p₁) ∝ S * vec(p₀)` (see [`advance!`](@ref)). ``S`` is stored as `transpose(Sᵀ)`, with a sparse or a dense `Sᵀ`
(see [`SparseMX`](@ref), [`DenseMX`](@ref)); `Matrix(S)` gives a dense copy.
"""
stepMX(PI::PathIntegration) = PI.stepMX[1]
stepMX(PI::PathIntegration, i) = PI.stepMX[i]

function PathIntegration(sde::AbstractSDE{d,k,m}, method::DiscreteTimeSteppingMethod, ts, axes::Vararg{Any,d}; kwargs...) where {d,k,m}
    sdestep = SDEStep(sde, method, ts; kwargs...)
    PathIntegration(sdestep,ts,axes...; kwargs...)
end
# PathIntegration{dynT, pdT, tsT, tpdMX_type, Tstp_idx, IKT, kwargT}
"""
    advance_till_converged!(PI; rtol = 1e-6, Tmax = nothing, check_dt = PI.ts[end] - PI.ts[1], atol = rtol*check_dt, check_iter = nothing, maxiter = 100_000)

Advance `PI` with [`advance!`](@ref) until the response PDF converges to the stationary PDF (or, for time-periodic systems, to the periodic PDF).
After every `check_dt` long interval the change ``ε = ∫|p(t) - p(t - check_dt)| dx`` is computed, and the iteration stops when ``ε ≤`` `atol`,
or after `Tmax` time (after `maxiter` steps if `Tmax = nothing`).

# Keyword Arguments
- `rtol = 1e-6`: tolerance of the change per unit time (`atol = rtol*check_dt`)
- `check_dt`: the time between the checks. The default is the time span of `PI.ts`: one time step for a time-invariant system,
  one period for a time-periodic system (so the PDFs are compared at the same phase).
- `check_iter = nothing`: the number of time steps between the checks (instead of `check_dt`)
- `Tmax = nothing`, `maxiter = 100_000`: the maximum time and the maximum number of steps

A constant time step is assumed. Returns `(PI, ε)`, where `ε` is the vector of the computed changes (its first element is a placeholder).
[`steady_state!`](@ref) computes the stationary PDF directly: it is usually faster, especially if the PDF converges slowly (a second eigenvalue of the step matrix near 1).
"""
function advance_till_converged!(PI::PathIntegration; rtol = 1e-6, Tmax = nothing, check_dt = PI.ts[end] - PI.ts[1], maxiter = 100_000, atol = rtol*check_dt, check_iter = nothing)
    _dt = PI.ts isa Number ? PI.ts : PI.ts[2] - PI.ts[1]
    if check_iter isa Nothing
        chk_itr = Int((check_dt + sqrt(eps(check_dt))) ÷_dt) - 1;
        # Assuming constant time step
    else
        chk_itr = check_iter - 1
    end
    if Tmax isa Nothing
        _maxiter = maxiter
    else
        _maxiter = Int((Tmax+sqrt(eps(Tmax))) ÷_dt);
        # Assuming constant time step
    end

    iter = zero(_maxiter)
    
    ϵ = [100*atol];

    # Compare the PDFs `check_dt` apart (a full period for time-periodic systems)
    p_prev = similar(PI.pdf.p)
    while ϵ[end] > atol && iter < _maxiter
        p_prev .= PI.pdf.p
        for _ in 1:chk_itr+1
            advance!(PI)
        end
        push!(ϵ,integrate_diff(PI.pdf,p_prev))
        iter = iter + chk_itr + 1
    end
    PI, ϵ
end


"""
    advance!(PI)

Advance the response PDF `PI.pdf` by one time step: ``p ← S p / ∫(S p)`` with the step matrix ``S`` of the next time interval of `PI.ts`
(the intervals are used cyclically), and increase `PI.t` by the length of the interval.
The marginal PDFs are not updated (see [`update_mPDFs!`](@ref)).
"""
function advance!(PI::PathIntegration)
    mass = _advance_to_temp!(PI.p_temp,PI)
    _corr_to_temp!(PI.pdf.p,PI.p_temp,mass)
    nothing
end
# p_temp = S p, returns ∫(S p)
function _advance_to_temp!(p_temp::tT,PI::PathIntegration{dynT}) where {tT<:AbstractArray{T,d},dynT<:AbstractSDEStep{d}} where {T,d}
    S = next_stepMX(PI)
    mass = dot(current_stepMX_wts(PI), vec(PI.pdf.p)) # ∫(S p) = dot(Sᵀw, p)
    mul!(vec(p_temp), S, vec(PI.pdf.p))
    PI.t = PI.t + get_PIdt(PI)
    mass
end
function _corr_to_temp!(res::tT,p_temp::tT,mass::Number) where {tT<:AbstractArray{T,d}} where {T,d}
    _I = 1/mass;
    @. res = p_temp * _I
    nothing
end

# Quadrature weights w of the grid: ∫p = dot(w, vec(p))
quadrature_weights(pdf::InterpolatedFunction) = [prod(w) for w in Iterators.product(map(ax -> ax.wts, pdf.axes)...)]
# Sᵀw for each step matrix S
stepMX_weights(::Nothing, pdf) = nothing
stepMX_weights(stepMX::AbstractVector{<:AbstractMatrix}, pdf) = [stepMX_weights(S, pdf) for S in stepMX]
stepMX_weights(S::AbstractMatrix{<:Number}, pdf) = transpose(S) * vec(quadrature_weights(pdf))

@inline current_stepMX_wts(PI::PathIntegration{dynT, pdT,tsT}) where {dynT, pdT,tsT<:Number} = PI.stepMX_wts
@inline current_stepMX_wts(PI::PathIntegration{dynT, pdT,tsT}) where {dynT, pdT,tsT<:AbstractArray} = PI.stepMX_wts[PI.step_idx]

@inline function next_stepMX(PI::PathIntegration{dynT, pdT,tsT}) where {dynT, pdT,tsT<:Number}
    PI.stepMX
end
@inline function next_stepMX(PI::PathIntegration{dynT, pdT,tsT}) where {dynT, pdT,tsT<:AbstractArray}
    PI.step_idx =  mod1(PI.step_idx + 1, length(PI.stepMX))
    PI.stepMX[PI.step_idx]
end

@inline function get_PIdt(PI::PathIntegration{dynT, pdT,tsT}) where {dynT, pdT,tsT<:Number}
    PI.ts
end
@inline function get_PIdt(PI::PathIntegration{dynT, pdT,tsT}) where {dynT, pdT,tsT<:AbstractArray}
    PI.ts[PI.step_idx+1]-PI.ts[PI.step_idx]
end
# Computations utilites
function init_DiagonalNormalPDF(axes...; μ_init = nothing, σ_init = nothing, kwargs...)
    if μ_init isa Nothing
        μs = [(axis[end]+axis[1])/2 for axis in axes]
    elseif μ_init isa Number
        μs = [μ_init for _ in axes]
    elseif μ_init isa Union{NTuple{length(axes),<:Number},AbstractVector{<:Number}}
        @assert length(μ_init) == length(axes) "Wrong number of initial μ values are given"
        μs = μ_init
    end

    if σ_init isa Nothing
        σ²s = [((axis[end]-axis[1])/12)^2 for axis in axes]
    elseif σ_init isa Number
        σ²s = [σ_init^2 for _ in axes]
    elseif σ_init isa Union{NTuple{length(axes),<:Number},AbstractVector{<:Number}}
        @assert length(σ_init) == length(axes) "Wrong number of initial σ values are given"
        σ²s = σ_init.^2
    end

    DiagonalNormalPDF(μs, σ²s)
end

(f::DiagonalNormalPDF)(x...) = prod(normal1D_σ2(μ, σ², _x) for (μ, σ², _x) in zip(f.μ, f.σ², x))

# (`Vararg{Any,N} where N`: arguments that are only passed on are otherwise not specialised on)
(PI::PathIntegration)(x::Vararg{Any,N}) where N = PI.pdf(x...)

## Recompute functions
"""
    reinit_PI_pdf!(PI, f = nothing; reset_t = true, reset_step_index = true)

Reinitialise the response PDF of `PI` with the values of the function `f(x_1, …, x_d)` at the grid nodes, or with the initial diagonal Gaussian
of [`PathIntegration`](@ref) (`μ_init`, `σ_init`) if `f = nothing`. `f` does not need to be normalised ([`advance!`](@ref) normalises the PDF).

- `reset_t = true`: set `PI.t = 0`
- `reset_step_index = true`: the next [`advance!`](@ref) uses the first step matrix
"""
function reinit_PI_pdf!(PI::PathIntegration,f = nothing; reset_t= true, reset_step_index = true)
    if f isa Nothing
        _f = init_DiagonalNormalPDF(PI.pdf.axes...; PI.IK.kwargs...)
    elseif f isa Function
        _f = f
    end
    recycle_interpolatedfunction!(PI.pdf, _f)

    if reset_t
        PI.t = zero(PI.t)
    end
    if reset_step_index
        PI.step_idx = zero(PI.step_idx)
    end
    PI
end

"""
    recompute_PI!(PI; par = nothing, t = nothing, f = nothing, Q_reinit_pdf = false, Q_recompute_stepMX = true, reset_t = true, reset_step_index = true, rowcomputation = nothing)

Reinitialise the response PDF (if `Q_reinit_pdf`, with [`reinit_PI_pdf!`](@ref)`(PI, f)`) and recompute the step matrices
(if `Q_recompute_stepMX`, with [`recompute_stepMX!`](@ref) and the keyword arguments `par`, `t` and `rowcomputation`).

# Example
Stationary PDFs of the Duffing oscillator of [`SDE`](@ref) (`p = [ζ, λ, σ]`) for several damping ratios:
```julia
PI = PathIntegration(sde, RK4(), 0.02, QuinticAxis(-4., 4., 41), QuinticAxis(-4., 4., 41))
pdfs = map(0.1:0.1:0.5) do ζ
    recompute_PI!(PI; par = [ζ, 0.25, 1.0], Q_reinit_pdf = true)
    steady_state!(PI)
    copy(PI.pdf.p)
end
```
"""
function recompute_PI!(PI::PathIntegration; par = nothing, t = nothing, f = nothing, Q_reinit_pdf = false, reset_t= true, reset_step_index = true, Q_recompute_stepMX = true, rowcomputation = nothing)
    if Q_reinit_pdf
        reinit_PI_pdf!(PI, f)
    end
    if Q_recompute_stepMX
        recompute_stepMX!(PI, par = par, t = t, reset_t = reset_t, reset_step_index = reset_step_index, rowcomputation = rowcomputation)
    end
end
"""
    recompute_stepMX!(PI; par = nothing, t = nothing, reset_t = true, reset_step_index = true, rowcomputation = nothing)

Recompute the step matrices of `PI` in place, e.g. for new parameters or time steps. The compiled time stepping and the buffers are reused,
so this is much faster than a new [`PathIntegration`](@ref). The response PDF is not changed (see [`reinit_PI_pdf!`](@ref), [`recompute_PI!`](@ref)).

# Keyword Arguments
- `par = nothing`: new parameter values, copied into the parameters of the SDE (which have to be a mutable container of the same length)
- `t = nothing`: new time points (a vector or a range: a step matrix for each interval) or a time step (a number: a step matrix for `[0, t]`)
- `reset_t = true`: set `PI.t = 0`
- `reset_step_index = true`: the next [`advance!`](@ref) uses the first step matrix
- `rowcomputation = nothing`: the row computation (e.g. [`ThreadedRowComputation`](@ref)); `nothing` uses the one given to [`PathIntegration`](@ref)
"""
function recompute_stepMX!(PI::PathIntegration; par = nothing, t = nothing, reset_t= true, reset_step_index = true, rowcomputation = nothing)
    if !(par isa Nothing)
        PI.IK.sdestep.sde.par .= par;
    end

    if t isa AbstractVector
        if length(PI.IK.t) != length(t)
            resize!(PI.IK.t,length(t))
        end
        PI.IK.t .= t
    elseif t isa Number
        resize!(PI.IK.t,2);
        PI.IK.t[1] = zero(eltype(PI.IK.t))
        PI.IK.t[2] = t
    end
    resize_stepMX!(PI.stepMX, length(PI.IK.t) - 1)

    reinit_stepMX!(PI.stepMX)
    kwargs = rowcomputation isa Nothing ? PI.IK.kwargs : merge(PI.IK.kwargs, (; rowcomputation = rowcomputation))
    fill_stepMX_ts!(PI.stepMX, PI.IK; kwargs...)
    PI.stepMX_wts = stepMX_weights(PI.stepMX, PI.pdf)

    if reset_t
        PI.t = zero(PI.t)
    end
    if reset_step_index
        PI.step_idx = zero(PI.step_idx)
    end
    nothing
end

function resize_stepMX!(stepMX::AbstractVector, n)
    l = length(stepMX)
    if n < l
        resize!(stepMX, n)
    elseif n > l
        append!(stepMX, [deepcopy(stepMX[1]) for _ in l+1:n])
    end
    stepMX
end

function reinit_stepMX!(stepMX::SparseArrays.AbstractSparseMatrixCSC) where T
    resize!(stepMX.nzval,0)
    resize!(stepMX.rowval,0)
    fill!(stepMX.colptr,one(eltype(stepMX.colptr)))
end
function reinit_stepMX!(stepMX::ThreadedSparseMatrixCSC) where T
    reinit_stepMX!(stepMX.A)
end
function reinit_stepMX!(stepMX::Transpose) where T
    reinit_stepMX!(stepMX.parent)
end
function reinit_stepMX!(stepMX::AbstractMatrix{T}) where T
    fill!(stepMX,zero(T))
end
function reinit_stepMX!(stepMX::AbstractVector{amT}) where amT<:AbstractMatrix{T} where T
    foreach(reinit_stepMX!, stepMX) # (fill! on a wrapped sparse matrix sets every element)
end


"""
    update_mPDFs!(PI; detached = false)

Compute the marginal PDFs `PI.marginal_pdfs` from the current response PDF `PI.pdf`, by integrating out the other coordinates with the quadrature weights of the axes.
The marginal PDFs are requested with the `mPDF_IDs` keyword argument of [`PathIntegration`](@ref); [`advance!`](@ref) does not update them.

`PI.marginal_pdfs` is a marginal PDF, or a `Tuple` of them for several `mPDF_IDs`. A marginal PDF `mpdf` is evaluated as `mpdf(x...)`,
and its [`InterpolatedFunction`](@ref) is `mpdf.pdf`. Use `detached = true` if `PI` was loaded from a file (e.g. with JLD2.jl).

# Example
```julia
PI = PathIntegration(sde, RK4(), 0.02, QuinticAxis(-4., 4., 41), QuinticAxis(-4., 4., 41); mPDF_IDs = [1, 2])
steady_state!(PI)
update_mPDFs!(PI)
p_x, p_v = PI.marginal_pdfs
p_x(0.5)
```
"""
update_mPDFs!(PI::PathIntegration{dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT,kwargT}; kwargs...) where {dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT<:Nothing,kwargT} = nothing

function update_mPDFs!(PI::PathIntegration{dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT,kwargT}; kwargs...) where {dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT<:MarginalPDF,kwargT}
    update_mPDF!(PI.marginal_pdfs,PI.pdf; kwargs...)
    nothing
end
function update_mPDFs!(PI::PathIntegration{dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT,kwargT}; kwargs...) where {dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT,kwargT}
    update_mPDF!.(PI.marginal_pdfs,Ref(PI.pdf); kwargs...)
    nothing
end