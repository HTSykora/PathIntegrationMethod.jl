"""
    PathIntegration(sde, method, ts, axes...; kwargs...)
    
Compute a `PathIntegration` object for computing the response probability density function evolution of the stochastic dynamical system defined by `sde`.

# Arguments
- `sde::AbstractSDE{d,k,m}`: ``d``-dimensional dynamic system which is subjected to an ``m``-dimensional Wiener process acting on the ``k...d`` coordinates.
- `method::DiscreteTimeSteppingMethod`: method used for the time-discretisation of the dynamical system described by `sde`. 
- `ts::Union{Number,AbstractArray{<:Number}}`: If `AbstractArray{<:Number}` is provided then the `stepMX`-s between the times in `ts` are computed. If `Number` is provided, then it computes a single `stepMX` is computed between 0 and `ts`.
- `axes::Vararg{aT,d} where aT<:GridAxis`: The axes spanning the space where the transitional PDF is computed for the dynamical system `sde`.

# Keyword Arguments
- `discreteintegrator = defaultdiscreteintegrator(sde, di_N = 31)`: Discrete integrator to evaluate the Chapman-Kolmogorov equation
    - The default discrete integrator algorithms are
        - `defaultdiscreteintegrator(sde::AbstractSDE{d,k,m}, di_N = 31) = GaussLegendreIntegrator(di_N)`
        - `defaultdiscreteintegrator(sde::SDE_VIO, di_N = 31) = Tuple(GaussLegendreIntegrator(di_N) for _ in 1:2)`
    - `di_N = 31`: Resolution of the discrete integrator. Can be a `Integer` or `NTuple{d-k+1,<:Integer}` that defines the discrete integration resolution in each `d-k+1` integration direction.
- `smart_integration = true`: Only integrate where the transitional PDF has nonzero elements. It is approximated with the step function. Use `false` if (time step * diffusion) results in a wide TPDF. Usually `true` is the better choice.
- `int_limit_thickness_multiplier = 6`: The "thickness" scaling of the TPDF during smart integration.
- `initialise_pdf = true`: Initialize the response probability density function. If false, then the RPDF is initialzed as p(x) ≡ 0.
- `f_init = nothing`: Initial RPDF as a function. If `f_init` is a `Nothing` and `initialise_pdf = true` then a diagonal Gaussian distribution is used as initial distribution
- `μ_init = nothing`: Mean of the initial Gaussian distribution used in case `f_init = nothing`. 
    - `Nothing`: Uses the middle of the range defined by the axis
    - `Number`: Uses the single value `μ_init` for each axis direction
    - `Union{NTuple{d,<:Number},AbstractVector{<:Number}}`: Individual means for each axis.
- `σ_init = nothing`: Standard deviation of the initial diagonal Gaussian distribution used in case `f_init = nothing`. 
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
function PathIntegration(sdestep::AbstractSDEStep{d,k,m}, _ts, axes::Vararg{Any,d}; 
    di_N = 31, discreteintegrator = defaultdiscreteintegrator(sdestep.sde, di_N),
    initialise_pdf = true, f_init = nothing, pre_compute = true, stepMXtype = nothing, sparse_tol = 1e-6,
    mPDF_IDs = nothing, extract_IK = Val{false}(), kwargs...) where {d,k,m}
    if stepMXtype isa StepMatrixRepresentation
        _stepMXtype = stepMXtype
    else
        _stepMXtype = get_stepMXtype(sdestep.sde, get_val_itp_type(axes); sparse_tol = sparse_tol, kwargs...)
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
        # * for CPU parallelisation extend this
        itpVs = Tuple(zero(axis.temp) for axis in axes); # = itpVs;
        
        ikt = IK_temp(BI_product(_eachindex.(itpVs)...), # = idx_it
                        BI_product(_val.(itpVs)...), # = val_it
                        itpVs, # = itpVs
                        zero(pdf.p)) # = itpM

        di = DiscreteIntegrator(discreteintegrator, sdestep, pdf.p, axes[k:end]...; kwargs...)
        IK = IntegrationKernel(sdestep, nothing, di, ts, pdf, ikt, (;sparse_tol = get_tol(_stepMXtype), kwargs...))
        
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
    PathIntegration(sdestep, pdf, p_temp,ts, stepMX, step_idx, IK, mpdf, kwargs,t)
end
_val(vals) = vals
get_ts(_ts::AbstractVector{tsT}) where tsT<:Number = _ts
get_ts(_ts::tsT) where tsT<:Number = [zero(_ts), _ts]
stepMX(PI::PathIntegration) = PI.stepMX[1]
stepMX(PI::PathIntegration, i) = PI.stepMX[i]

function PathIntegration(sde::AbstractSDE{d,k,m}, method::DiscreteTimeSteppingMethod, ts, axes::Vararg{Any,d}; kwargs...) where {d,k,m}
    sdestep = SDEStep(sde, method, ts; kwargs...)
    PathIntegration(sdestep,ts,axes...; kwargs...)
end
# PathIntegration{dynT, pdT, tsT, tpdMX_type, Tstp_idx, IKT, kwargT}
function advance_till_converged!(PI::PathIntegration; rtol = 1e-6, Tmax = nothing, check_dt = PI.ts[end], maxiter = 100_000, atol = rtol*check_dt, check_iter = nothing)
    _dt = PI.ts isa Number ? PI.ts : PI.ts[2]
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

    while ϵ[end] > atol && iter < _maxiter
        for _ in 1:chk_itr
            advance!(PI)
        end
        _advance_to_temp!(PI.p_temp,PI)
        _corr_to_temp!(PI.p_temp,PI.p_temp,PI)
        push!(ϵ,integrate_diff(PI.pdf,PI.p_temp))
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
function _advance_to_temp!(p_temp::tT,PI::PathIntegration{dynT}) where {tT<:AbstractArray{T,d},dynT<:AbstractSDEStep{d}} where {T,d}
    mul!(vec(p_temp), next_stepMX(PI), vec(PI.pdf.p))
    PI.t = PI.t + get_PIdt(PI)
    nothing
end
function _corr_to_temp!(res::tT,p_temp::tT,PI::PathIntegration{dynT}) where {tT<:AbstractArray{T,d},dynT<:AbstractSDEStep{d}} where {T,d}
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

(PI::PathIntegration)(x...) = PI.pdf(x...)

## Recompute functions
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

function recompute_PI!(PI::PathIntegration; par = nothing, t = nothing, f = nothing, Q_reinit_pdf = false, reset_t= true, reset_step_index = true, Q_recompute_stepMX = true)
    if Q_reinit_pdf
        reinit_PI_pdf!(PI, f)
    end
    if Q_recompute_stepMX
        recompute_stepMX!(PI, par = par, t = t, reset_t = reset_t, reset_step_index = reset_step_index)
    end
end
function recompute_stepMX!(PI::PathIntegration; par = nothing, t = nothing, reset_t= true, reset_step_index = true)
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
    fill!.(stepMX,zero(T))
end


update_mPDFs!(PI::PathIntegration{dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT,kwargT}; kwargs...) where {dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT<:Nothing,kwargT} = nothing

function update_mPDFs!(PI::PathIntegration{dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT,kwargT}; kwargs...) where {dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT<:MarginalPDF,kwargT}
    update_mPDF!(PI.marginal_pdfs,PI.pdf; kwargs...)
    nothing
end
function update_mPDFs!(PI::PathIntegration{dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT,kwargT}; kwargs...) where {dynT, pdT, tsT, stepmxT, Tstp_idx, IKT, ptempT,mpdtT,kwargT}
    update_mPDF!.(PI.marginal_pdfs,Ref(PI.pdf); kwargs...)
    nothing
end