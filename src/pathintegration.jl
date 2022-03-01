"""
    PathIntegration(sde, method, ts, axes...; kwargs...)
    
Compute a `PathIntegration` object for computing the response probability density function evolution of the stochastic dynamical system defined by `sde`.

# Arguments
- `sde::AbstractSDE{d,k,m}`: ``d``-dimensional dynamic system which is subjected to an ``m``-dimensional Wiener process acting on the ``k...d`` coordinates.
- `method::DiscreteTimeSteppingMethod`: method used for the time-discretisation of the dynamical system described by `sde`. 
- `ts::Union{Number,AbstractArray{<:Number}}`: If `AbstractArray{<:Number}` is provided then the `stepMX`-s between the times in `ts` are computed. If `Number` is provided, then it computes a single `stepMX` is computed between 0 and `ts`.
- `axes::Vararg{aT,d} where aT<:GridAxis`: The axes spanning the space where the transitional PDF is computed for the dynamical system `sde`.

# Keyword Arguments
- `discreteintegrator = d==k ? QuadGKIntegrator() : ClenshawCurtisIntegrator()`: Discrete integrator to evaluate the Chapman-Kolmogorov equation
- `di_N = 21`: Resolution of the discrete integrator. Can be an `NTuple{d-k+1,<:Integer}` or an `AbstractArray{<:Integer}` to define resolution to each integration variable
- `smart_integration = true`: Only integrate where the transitional PDF has nonzero elements. It is approximated with the step function. Use `false` if (time step * diffusion) results in a wide TPDF. Usually `true` is the better choice.
- `int_limit_thickness_multiplier = 6`: The "thickness" scaling of the TPDF during smart integration.
- `initialise_pdf = true`: Initialize the response probability density function. If false, then the RPDF is initialzed as p(x) ≡ 0.
- `f_init = nothing`: Initial RPDF as a function. If `f_init` is a `Nothing` and `initialise_pdf = true` then a diagonal Gaussian distribution is used as initial distribution
- `μ_init = nothing`: Mean of the initial Gaussian distribution. 
    - `Nothing`: Uses the middle of the range defined by the axis
    - `Number`: Uses the single value `μ_init` for each axis direction
    - `Union{NTuple{d,<:Number},AbstractVector{<:Number}}`: Individual means for each axis.
- `σ_init = nothing`: Standard deviation of the initial diagonal Gaussian distribution. 
    - `Nothing`: Uses the 1/12th of the range width in each direction defined by the `axes`.
    - `Number`: Uses the single value `σ_init` for each axis direction
    - `Union{NTuple{d,<:Number},AbstractVector{<:Number}}`: Individual standard deviations for each axis.
- `pre_compute = true`: Compute the `stepMX`. This should be left unchanged if the RPDF computation is the goal.
- `stepMXtype = nothing`: Step matrix representation type. The default representation depends on `d` and the interpolation used (see Interpolation): `d ≤ 2` with sparse interpolation the `stepMX` is a multithreaded sparse matrix, for dense interpolations it is dense, and for `d>2` the default is a multithreaded sparse
    Possible options:
        - `SparseMX(; threaded = true, sparse_tol = 1e-6)`
        - `DenseMX()`
- `multithreaded_sparse = true`: is the sparse `stepMX` multithreaded if no `stepMXtype` is specified
- `sparse_tol = 1e-6`: absolute tolerance for the elements considered as zero values in the sparse stepMX if no `stepMXtype` is specified
- `mPDF_IDs = nothing`: Marginal PDF (mPDF) for IDinates specified by `mPDF_IDs`
    - `Nothing`: No mPDF is initialised.
    - `Integer`: 1-dimensional mPDF is initalised for state-variable `mPDF_IDs`
    - `NTuple{n,<:Integer}`: `n`-dimensional mPDF is initalised for state-variables specified in `mPDF_IDs`
    - `Union{Tuple{NTuple{n,<:Integer}},Vector{NTuple{d,<:Integer}}}`: Multiple mPDF initialised
- `allow_extrapolation::Bool = false`: Allow nonzero extrapolation outside the region specified by `axes`
- `zero_extrapolation::Bool = true`: Return 0 if extrapolated outside of the region specified by `axes`
----
For methods, discrete integrators, interpolators, and examples please refer to the documentation. 
"""
function PathIntegration(sde::AbstractSDE{d,k,m}, method, ts, axes::Vararg{Any,d}; 
    discreteintegrator = d==k ? QuadGKIntegrator() : ClenshawCurtisIntegrator(),
    di_N = 21,  # discrete integration resolution
    initialise_pdf = true, f_init = nothing, pre_compute = true, stepMXtype = nothing, sparse_tol = 1e-6,
    mPDF_IDs = nothing, extract_IK = Val{false}(), kwargs...) where {d,k,m}
    if stepMXtype isa StepMatrixRepresentation
        _stepMXtype = stepMXtype
    else
        _stepMXtype = get_stepMXtype(sde, get_val_itp_type(axes); sparse_tol = sparse_tol, kwargs...)
    end

    if method isa DiscreteTimeSteppingMethod
        x0 = zeros(d) # * not type safe for autodiff
        x1 = similar(x0)
        sdestep = SDEStep(sde, method, x0, x1, ts; kwargs...)
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
    pdf = InterpolatedFunction(axes...; f = _f, kwargs...)

    step_idx = 0
    t = 0.
    if pre_compute
        # * for CPU parallelisation extend this
        itpVs = Tuple(zero(axis.temp) for axis in axes); # = itpVs;
        
        ikt = IK_temp(BI_product(_eachindex.(itpVs)...), # = idx_it
                        BI_product(_val.(itpVs)...), # = val_it
                        itpVs, # = itpVs
                        zero(pdf.p)) # = itpM
        if di_N isa Number
            di_res = Tuple(di_N for _ in 1:k-d+1)
        elseif di_N isa NTuple{d-k+1,<:Integer}
            di_res = di_N
        elseif di_N isa AbstractArray{<:Integer}
            if length(di_res) == d-k+1
                di_res = di_N
            end
        else
            error("di_N is not a Union{Number, NTuple{M,Integer}, AbstractArray{Integer}} with length $(d-k+1)")
        end
        di = DiscreteIntegrator(discreteintegrator, pdf.p, di_res, axes[k:end]...; kwargs...)
        IK = IntegrationKernel(sdestep, nothing, di, ts, pdf, ikt, kwargs)
        
        if extract_IK isa Val{true}
            return IK
        end
        stepMX = compute_stepMX(IK; stepMXtype = _stepMXtype, kwargs...)
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
    PathIntegration(sdestep, pdf, p_temp,ts, stepMX, step_idx, IK, mpdf, kwargs,t)
end
_val(vals) = vals
# PathIntegration{dynT, pdT, tsT, tpdMX_type, Tstp_idx, IKT, kwargT}
function advance_till_converged!(PI::PathIntegration; rtol = 1e-6, Tmax = nothing, check_dt = PI.ts[1], check_iter = nothing , maxiter = 100_000, atol = rtol*check_dt)
    if check_iter isa Nothing
        chk_itr = Int((check_dt+sqrt(eps(check_dt))) ÷ PI.ts[1]) - 1;
        # Assuming constant time step
    else
        chk_itr = check_iter - 1
    end
    if Tmax isa Nothing
        _maxiter = maxiter
    else
        _maxiter = Int((Tmax+sqrt(eps(Tmax))) ÷ PI.ts[1]);
        # Assuming constant time step
    end

    iter = zero(_maxiter)
    
    ϵ = 100*atol;

    while ϵ > atol && iter < _maxiter
        for _ in 1:chk_itr
            advance!(PI)
        end
        _advance_to_temp!(PI.p_temp,PI)
        _corr_to_temp!(PI.p_temp,PI.p_temp,PI)
        ϵ = integrate_diff(PI.pdf,PI.p_temp)
        @. PI.pdf.p = PI.p_temp
        iter = iter + chk_itr + 1
    end
    PI, ϵ
end


function advance!(PI::PathIntegration)
    _advance_to_temp!(PI.p_temp,PI)
    _corr_to_temp!(PI.pdf.p,PI.p_temp,PI)
    nothing
end
function _advance_to_temp!(p_temp::tT,PI::PathIntegration{dynT}) where {tT<:AbstractArray{T,d},dynT<:SDEStep{d}} where {T,d}
    mul!(vec(p_temp), next_stepMX(PI), vec(PI.pdf.p))
    PI.t = PI.t + get_PIdt(PI)
    nothing
end
function _corr_to_temp!(res::tT,p_temp::tT,PI::PathIntegration{dynT}) where {tT<:AbstractArray{T,d},dynT<:SDEStep{d}} where {T,d}
    _I = 1/_integrate(p_temp, PI.pdf.axes...);
    @. res = PI.p_temp * _I
    nothing
end

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
    PI.ts[PI.step_idx+1]
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

(PI::PathIntegration)(x...) = PI.pdf(x...)

## Recompute functions
function reinit_PI_pdf!(PI::PathIntegration,f = nothing, reset_t= true, reset_step_index = true)
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
end


function recompute_stepMX!(PI::PathIntegration; par = nothing, t = nothing, f = nothing, Q_reinit_pdf = false, reset_t= true, reset_step_index = true)
    if !(par isa Nothing)
        PI.IK.sdestep.sde.par .= par;
    end

    if t isa AbstractVector
        if length(PI.IK.ts) != length(ts)
            resize!(PI.IK.ts,length(ts))
        end
        PI.IK.t .= ts
    elseif t isa Number
        resize!(PI.IK.t,2);
        PI.IK.t[1] = zero(eltype(PI.IK.t))
        PI.IK.t[2] = t
    end
    if Q_reinit_pdf
        reinit_PI_pdf!(PI, f)
    end

    reinit_stepMX!(PI.stepMX)
    fill_stepMX_ts!(PI.stepMX, PI.IK; PI.IK.kwargs...)

    if reset_t
        PI.t = zero(PI.t)
    end
    if reset_step_index
        PI.step_idx = zero(PI.step_idx)
    end
    nothing
end

function reinit_stepMX!(stepMX::AbstractMatrix{T}) where T
    fill!(stepMX,zero(T))
end
function reinit_stepMX!(stepMX::AbstractVector{amT}) where amT<:AbstractMatrix{T} where T
    fill!.(stepMX,zero(T))
end


update_mPDFs!(PI::PathIntegration{dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT,kwargT}) where {dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT<:Nothing,kwargT} = nothing

function update_mPDFs!(PI::PathIntegration{dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT,kwargT}) where {dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT<:MarginalPDF,kwargT}
    update_mPDF!(PI.marginal_pdfs,PI.pdf)
    nothing
end
function update_mPDFs!(PI::PathIntegration{dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT,kwargT}) where {dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT,kwargT}
    update_mPDF!.(PI.marginal_pdfs,Ref(PI.pdf))
    nothing
end